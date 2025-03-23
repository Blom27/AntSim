using UnityEngine;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;
using Unity.Burst;
using System.Collections.Generic;
using Unity.Collections.LowLevel.Unsafe;
using System.Linq;

public class AntSimulationController : MonoBehaviour {
    [Header("Optimisation Settings")]
    public bool drawPheromones = true;
    [Range(1, 2)]
    public float simulationSpeed = 1;
    public int pheromoneThreshold = 10000;
    [Range(float.Epsilon, 1 - float.Epsilon)]
    public float pheromoneTrimPercentage = 0.1f;
    
    [Header("Ant Pathfinding Settings")]
    [Range(0.0001f, 1)]
    public float pheromoneAgeSensitivity;
    [Range(0.0001f, 1)]
    public float pheromoneTypeSensitivity;
    [Range(0.0001f, 1)]
    public float distanceHomeSensitivity;
    [Range(0.0001f, 1)]
    public float directionSensitivity;
    [Range(0.0001f, 1)]
    public float sensorDistanceSensitivity;
    
    [Header("Simulation Settings")]
    public int colonySize = 50;
    public float2 nestPosition = Vector2.zero;
    public float worldSize = 50f;
    public RenderMode currentRenderMode = RenderMode.Gizmos;
    public ISimulationRenderer _activeRenderer;

    public enum RenderMode { Gizmos/*, Sprites, ThreeD*/ }
    
    [Header("Data")]
    public AntSimulationData _simulationData;

    private NativeArray<PheromoneData> _pheromoneJobData;

    void Start() {
        pheromoneAgeSensitivity = 1;
        pheromoneTypeSensitivity = 1;
        distanceHomeSensitivity = 1;
        directionSensitivity = 1;
        sensorDistanceSensitivity = 1;
        _simulationData = AntSimulationCore.Initialize(_simulationData, nestPosition, 3, worldSize, 0.5f, new MovementParamaters {
            pheromoneAgeSensitivity = pheromoneAgeSensitivity,
            pheromoneTypeSensitivity = pheromoneTypeSensitivity,
            distanceHomeSensitivity = distanceHomeSensitivity,
            directionSensitivity = directionSensitivity,
            sensorDistanceSensitivity = sensorDistanceSensitivity
        });
        _simulationData.PheromoneThreshold = pheromoneThreshold;
        if(colonySize > 256) CreateAntsWithJobSystem(colonySize);
        else _simulationData = AntSimulationCore.CreateAnts(colonySize, _simulationData);
        SwitchRenderer(currentRenderMode);
    }

    void Update() {
        if (_simulationData.Ants.Count > 0) {
            _simulationData.movementSettings = new MovementParamaters {
                pheromoneAgeSensitivity = pheromoneAgeSensitivity,
                pheromoneTypeSensitivity = pheromoneTypeSensitivity,
                distanceHomeSensitivity = distanceHomeSensitivity,
                directionSensitivity = directionSensitivity,
                sensorDistanceSensitivity = sensorDistanceSensitivity
            };
            _simulationData.PheromoneThreshold = pheromoneThreshold;
            if (Input.GetMouseButtonDown(0)) _simulationData.FoodSources.Add(AntSimulationCore.CreateNewFoodSource(_simulationData));
            if(_simulationData.Pheromones.Count > 64){
                _simulationData = AntSimulationCore.LocationCheck(_simulationData);
                _simulationData.SimulationTime = Time.deltaTime * simulationSpeed;
                _simulationData = AntSimulationCore.UpdateFoodSources(_simulationData);
                _simulationData = AntSimulationCore.UpdateAntEnergy(_simulationData);
                UpdatePheromonesWithJobSystem(Time.deltaTime * simulationSpeed);
                _simulationData = AntSimulationCore.UpdatePheromoneTimers(_simulationData);
                _simulationData = AntSimulationCore.UpdateAntDirections(_simulationData);
                _simulationData = AntSimulationCore.ClearPheromonesInRangeAndPositionSensors(_simulationData);
                UpdatePheromonesInRangeWithJobSystem();
                _simulationData = AntSimulationCore.UpdateAntRotations(_simulationData);
                _simulationData = AntSimulationCore.MoveAnts(_simulationData); // Move ants
            }
            else _simulationData = AntSimulationCore.Update(_simulationData, Time.deltaTime * simulationSpeed);
            _simulationData = AntSimulationCore.OrganiseLists(_simulationData);
        }
        _simulationData = AntSimulationCore.AntCountCheck(_simulationData);
        if(_simulationData.Pheromones.Count >= pheromoneThreshold * 1.1) {
            for(int i = 0; i < pheromoneThreshold * pheromoneTrimPercentage; i++) {
                _simulationData.Pheromones.RemoveAt(0);
            } 
        }
        if(currentRenderMode != RenderMode.Gizmos) _activeRenderer.Draw(_simulationData, drawPheromones);
    }

    void OnDrawGizmos() {
        if(currentRenderMode == RenderMode.Gizmos && Application.isPlaying) {
            _activeRenderer.Draw(_simulationData, drawPheromones);
        }
    }
    
    void OnDestroy() {
        if (_pheromoneJobData.IsCreated) _pheromoneJobData.Dispose();
    }

    public void SwitchRenderer(RenderMode newMode) {
        _activeRenderer?.CleanUp();

        switch(newMode) {
            case RenderMode.Gizmos:
                _activeRenderer = new GizmoRenderer();
                break;
                /*
            case RenderMode.Sprites:
                //_activeRenderer = new SpriteRenderer();
                break;
            case RenderMode.ThreeD:
                // _activeRenderer = new ThreeDRenderer();
                break;
                */
        }
        _activeRenderer.Initialize();
        currentRenderMode = newMode;
    }

    public void CreateAntsWithJobSystem(int colonySize) {
        // Create a single native array for all ant data
        NativeArray<AntJobData> antJobData = new NativeArray<AntJobData>(colonySize, Allocator.TempJob);
        
        // Set up and schedule the job
        var createAntsJob = new CreateAntsJob {
            NestPosition = new float2(nestPosition.x, nestPosition.y),
            NestSize = _simulationData.NestSize,
            Ants = antJobData
        };
        
        // Execute the job
        JobHandle jobHandle = createAntsJob.Schedule(colonySize, 64);
        jobHandle.Complete();
        
        // Transfer results to AntSimulationCore
        _simulationData.Ants = antJobData.Select(ajd => new AntData {
            baseData = ajd,
            pheromoneIndiciesInRange = new HashSet<int>(),
            foodCollected = 0
        }).ToList();
        _simulationData.ColonySize = _simulationData.Ants.Count;

        antJobData.Dispose();
    }

    private void UpdatePheromonesWithJobSystem(float deltaTime) {
        int pheromoneCount = _simulationData.Pheromones.Count;
        if (pheromoneCount == 0) return;
        // Initialize native containers
        using (NativeQueue<int> removeQueue = new NativeQueue<int>(Allocator.TempJob)) {
            // Resize native array if needed
            if (!_pheromoneJobData.IsCreated || _pheromoneJobData.Length < pheromoneCount) {
                if (_pheromoneJobData.IsCreated) _pheromoneJobData.Dispose();
                _pheromoneJobData = new NativeArray<PheromoneData>(pheromoneCount, Allocator.Persistent);
            }
            // Copy data to native array
            for (int i = 0; i < pheromoneCount; i++) {
                _pheromoneJobData[i] = new PheromoneData {
                    position = new float2(_simulationData.Pheromones[i].position.x, 
                                        _simulationData.Pheromones[i].position.y),
                    owner = _simulationData.Pheromones[i].owner,
                    age = _simulationData.Pheromones[i].age
                };
            }
            var updateJob = new UpdatePheromoneAgeJob {
                PheromoneData = _pheromoneJobData,
                DeltaTime = deltaTime,
                RemoveQueueWriter = removeQueue.AsParallelWriter()
            };
            JobHandle jobHandle = updateJob.Schedule(pheromoneCount, 64);
            jobHandle.Complete();

            // Process expired pheromones
            NativeArray<int> indicesToRemove = removeQueue.ToArray(Allocator.Temp);// Create a HashSet of indices to remove for quick lookup
            HashSet<int> removeIndices = new HashSet<int>();
            for (int i = 0; i < indicesToRemove.Length; i++) removeIndices.Add(indicesToRemove[i]);

            // Create a new pheromone list from the updated job data
            List<PheromoneData> newPheromones = new List<PheromoneData>();
            for (int i = 0; i < pheromoneCount; i++) {
                // Skip if this pheromone should be removed
                if (removeIndices.Contains(i)) continue;
                
                // Add to new list with updated age value
                newPheromones.Add(new PheromoneData {
                    position = new Vector2(_pheromoneJobData[i].position.x, _pheromoneJobData[i].position.y),
                    owner = _pheromoneJobData[i].owner,
                    age = _pheromoneJobData[i].age
                });
            }
            // Replace the old list with the new one
            _simulationData.Pheromones = newPheromones;
            indicesToRemove.Dispose();
        }
    }
    
    private void UpdatePheromonesInRangeWithJobSystem() {
        int antCount = _simulationData.Ants.Count;
        int pheromoneCount = _simulationData.Pheromones.Count;
        if (antCount == 0 || pheromoneCount == 0) return;

        // Initialize native containers
        NativeArray<AntJobData> antJobData = new NativeArray<AntJobData>(antCount, Allocator.TempJob);
        NativeArray<PheromoneData> pheromoneJobData = new NativeArray<PheromoneData>(pheromoneCount, Allocator.TempJob);
        NativeArray<float2> sensorPositions = new NativeArray<float2>(antCount * 3, Allocator.TempJob); // 3 sensors per ant

        // Copy data to native arrays
        for (int i = 0; i < antCount; i++) {
            antJobData[i] = _simulationData.Ants[i].baseData;
            float2[] sensorPos = AntSimulationCore.CalculateSensorPositions(_simulationData.Ants[i].baseData);
            for (int j = 0; j < 3; j++) {
                sensorPositions[i * 3 + j] = sensorPos[j];
            }
        }

        for (int i = 0; i < pheromoneCount; i++) {
            pheromoneJobData[i] = _simulationData.Pheromones[i];
        }

        // Define chunk size (e.g., 64 pheromones per job)
        int chunkSize = 64;
        int numChunks = (int)(pheromoneCount / chunkSize) + 1;

        // Create arrays to store job handles and results
        NativeArray<JobHandle> jobHandles = new NativeArray<JobHandle>(numChunks, Allocator.Temp);
        NativeArray<UnsafeParallelMultiHashMap<int, int>> pheromoneMaps = new NativeArray<UnsafeParallelMultiHashMap<int, int>>(numChunks, Allocator.Temp);

        // Launch jobs for each chunk
        for (int chunkIndex = 0; chunkIndex < numChunks; chunkIndex++) {
            int chunkStartIndex = chunkIndex * chunkSize;
            int chunkEndIndex = math.min(chunkStartIndex + chunkSize, pheromoneCount);

            // Create a new map for this job
            pheromoneMaps[chunkIndex] = new UnsafeParallelMultiHashMap<int, int>(antCount, Allocator.TempJob);

            var processJob = new ProcessPheromonesJob {
                Ants = antJobData,
                Pheromones = pheromoneJobData,
                SensorPositions = sensorPositions,
                PheromonesInRangeMap = pheromoneMaps[chunkIndex],
                StartIndex = chunkStartIndex,
                EndIndex = chunkEndIndex
            };

            // Schedule the job
            jobHandles[chunkIndex] = processJob.Schedule();
        }

        // Wait for all jobs to complete
        JobHandle.CompleteAll(jobHandles);

        // Combine results from all jobs
        HashSet<int> pheromonesToRemove = new HashSet<int>();
        for (int chunkIndex = 0; chunkIndex < numChunks; chunkIndex++) {
            foreach (var pair in pheromoneMaps[chunkIndex]) {
                int antIndex = pair.Key;
                int pheromoneIndex = pair.Value;

                if (antIndex == -1) {
                    // This pheromone is marked for removal
                    if(_simulationData.Pheromones.Count >= pheromoneThreshold) pheromonesToRemove.Add(pheromoneIndex);
                } else {
                    // Add the pheromone to the ant's pheromonesInRange list
                    if (!_simulationData.Ants[antIndex].pheromoneIndiciesInRange.Contains(pheromoneIndex)) {
                        _simulationData.Ants[antIndex].pheromoneIndiciesInRange.Add(pheromoneIndex);
                    }
                }
            }
        }

        // Remove irrelevant pheromones
        List<PheromoneData> newPheromones = new List<PheromoneData>();
        for (int i = 0; i < pheromoneCount; i++) {
            if (!pheromonesToRemove.Contains(i)) {
                newPheromones.Add(_simulationData.Pheromones[i]);
            }
        }
        _simulationData.Pheromones = newPheromones;

        // Dispose native containers
        antJobData.Dispose();
        pheromoneJobData.Dispose();
        sensorPositions.Dispose();
        for (int i = 0; i < numChunks; i++) {
            pheromoneMaps[i].Dispose();
        }
        jobHandles.Dispose();
        pheromoneMaps.Dispose();
    }
}

[System.Serializable]
public class GizmoRenderer : ISimulationRenderer {
    private Color _antColor = Color.white;
    private Color _foodColor = Color.green;

    public void Initialize() { /* Not needed for Gizmos */ }

    public void Draw(AntSimulationData data, bool drawPheromones = true) {
        if (drawPheromones)DrawPheromones(data.Pheromones);
        DrawNest(data.NestPosition);
        DrawAnts(data.Ants);
        DrawFood(data.FoodSources);
    }

    private void DrawNest(Vector2 position) {
        Gizmos.color = Color.red;
        Gizmos.DrawWireSphere(position, 1f);
    }

    private void DrawAnts(List<AntData> ants) {
        foreach (var ant in ants) {
            Gizmos.color = (ant.baseData.isCarryingFood == 1) ? Color.yellow : _antColor;
            Gizmos.DrawCube(new Vector3(ant.baseData.position.x, ant.baseData.position.y, 0), Vector3.one * 0.2f);
        }
    }

    private void DrawPheromones(List<PheromoneData> pheromones) {
        foreach (var p in pheromones) {
            float ageLerp = Mathf.Clamp01(p.age / p.owner.maxEnergy);
            Gizmos.color = Color.Lerp(
                (p.owner.isCarryingFood == 1) ? Color.green : Color.blue, 
                Color.red, 
                ageLerp
            );
            Gizmos.DrawSphere(new Vector3(p.position.x, p.position.y, 0), 0.1f * (1 - ageLerp)); // Shrink with age
        }
    }

    private void DrawFood(List<FoodSource> foodSources) {
        Gizmos.color = _foodColor;
        foreach(var food in foodSources) {
            Gizmos.DrawWireSphere(food.position, (float)food.size);
        }
    }

    public void CleanUp() { /* Gizmos auto-clean */ }
}

public static class AntSimulationCore {

    public static Unity.Mathematics.Random random = new Unity.Mathematics.Random((uint)System.DateTime.Now.Millisecond);

    public static AntSimulationData Initialize(AntSimulationData simData, Vector2 nestPos, int foodSourceCount, float WorldSize, float NestSize, MovementParamaters MovementSettings) {
        simData = new AntSimulationData {
            SimulationTime = 0f,
            ColonySize = 0,
            NestPosition = nestPos,
            NestSize = NestSize,
            WorldSize = WorldSize,
            Ants = new List<AntData>(),
            Pheromones = new List<PheromoneData>(),
            FoodSources = new List<FoodSource>(),
            movementSettings = MovementSettings
        };
        // Create food sources
        for(int i = 0; i < foodSourceCount; i++) simData.FoodSources.Add(CreateNewFoodSource(simData));
        return simData;
    }

    public static AntSimulationData CreateAnts(int colonySize, AntSimulationData simData) {
        simData.ColonySize += colonySize;
        for (int i = 0; i < simData.ColonySize; i++) {
            AntJobData ant = CreateAnt(simData.NestPosition, simData.NestSize);
            simData.Ants.Add(new AntData {
                baseData = ant,
                pheromoneIndiciesInRange = new HashSet<int>(),
                sensorPositions = new float2[3],
                sensorStrengths = new float[3]
            });
        }
        return simData;
    }

    public static AntJobData CreateAnt(float2 nestPosition, float nestSize) {
        float2 randomOffset = random.NextFloat2(-nestSize / 2, nestSize / 2);
        float2 position = nestPosition + randomOffset;

        float rotation = random.NextFloat(0, 360f);
        float Speed = random.NextFloat(0.01f, 1) * 2;
        float RotationSpeed = random.NextFloat(90, 360);
        float SensorAngle = random.NextFloat(15, 90); // Randomize sensor angle
        float SensorRadius = random.NextFloat(0.5f, 1) * 2;
        float PheromoneIntervalValue = (random.NextFloat(0.05f, 1) * Speed);
        float PheromoneDecayRate = random.NextFloat(0.001f, 0.03f);
        float PheromoneMergeRadius = (random.NextFloat(0.5f, 1f) * 0.45f) * Speed;
        float MaxEnergy = ((((PheromoneDecayRate / 10) + PheromoneMergeRadius + Speed + (RotationSpeed / 360)) / (PheromoneIntervalValue * PheromoneDecayRate)) / (SensorAngle / 91)) * random.NextInt(1, 50);
        float EnergyConsumptionRate = ((MaxEnergy / 10) / ((PheromoneDecayRate * PheromoneIntervalValue) * Speed * (RotationSpeed / 360))) / 1000;

        return new AntJobData {
            position = position,
            rotation = rotation,
            isCarryingFood = 0,

            pheromoneTimer = 0f,
            pheromoneIntervalValue = PheromoneIntervalValue,
            pheromoneDecayRate = PheromoneDecayRate,
            pheromoneMergeRadius = PheromoneMergeRadius,

            sensorRadius = SensorRadius,
            sensorAngle = SensorAngle,
            smellRadius = (2 / SensorRadius),

            speed = Speed,
            rotationSpeed = RotationSpeed,

            maxEnergy = MaxEnergy,
            energy = MaxEnergy,
            energyConsumptionRate = EnergyConsumptionRate,
            pheromoneMaxAge = (PheromoneDecayRate * PheromoneIntervalValue) * (MaxEnergy / EnergyConsumptionRate)
        };
    }

    public static AntJobData CreateOffspring(AntSimulationData simData, int parentCount) {
        float2 randomOffset = random.NextFloat2(-simData.NestSize / 2, simData.NestSize / 2);
        float2 position = simData.NestPosition + randomOffset;
        float rotation = 0;
        float PheromoneIntervalValue = 0;
        float PheromoneDecayRate = 0;
        float PheromoneMergeRadius = 0;
        float Speed = 0;
        float RotationSpeed = 0;
        float SensorAngle = 0;
        float SensorRadius= 0;
        float MaxEnergy = 0;
        for(int i = 0; i < parentCount; i++) {
            int randomIndex = (random.NextFloat(0, 1) < 0.02f) ? random.NextInt(0, simData.Ants.Count / 10) : 0;
            rotation += simData.Ants[i + randomIndex].baseData.rotation;
            PheromoneIntervalValue += simData.Ants[i + randomIndex].baseData.pheromoneIntervalValue;
            PheromoneDecayRate += simData.Ants[i + randomIndex].baseData.pheromoneDecayRate;
            PheromoneMergeRadius += simData.Ants[i + randomIndex].baseData.pheromoneMergeRadius;
            Speed += simData.Ants[i].baseData.speed;
            RotationSpeed += simData.Ants[i + randomIndex].baseData.rotationSpeed;
            SensorAngle += simData.Ants[i + randomIndex].baseData.sensorAngle;
            SensorRadius += simData.Ants[i + randomIndex].baseData.sensorRadius;
            MaxEnergy += simData.Ants[i + randomIndex].baseData.maxEnergy;
        }
        rotation /= parentCount;
        PheromoneIntervalValue /= parentCount;
        PheromoneDecayRate /= parentCount;
        PheromoneMergeRadius /= parentCount;
        Speed /= parentCount;
        RotationSpeed /= parentCount;
        SensorAngle /= parentCount;
        SensorRadius /= parentCount;
        MaxEnergy /= parentCount;
        PheromoneMergeRadius *= SensorRadius;
        float SmellRadius = (2 / SensorRadius);
        float EnergyConsumptionRate = ((MaxEnergy / 10) / ((PheromoneDecayRate * PheromoneIntervalValue) * Speed * (RotationSpeed / 180))) / 1000;

        return new AntJobData {
            position = position,
            rotation = NormaliseAngle(rotation + 180),
            isCarryingFood = 0,

            pheromoneTimer = 0f,
            pheromoneIntervalValue = PheromoneIntervalValue,
            pheromoneDecayRate = PheromoneDecayRate,
            pheromoneMergeRadius = PheromoneMergeRadius,

            sensorRadius = SensorRadius,
            sensorAngle = SensorAngle,
            smellRadius = SmellRadius,

            speed = Speed,
            rotationSpeed = RotationSpeed,
            maxEnergy = MaxEnergy,
            energy = MaxEnergy,
            energyConsumptionRate = EnergyConsumptionRate,
            pheromoneMaxAge = (PheromoneDecayRate * PheromoneIntervalValue) * (MaxEnergy / EnergyConsumptionRate)
        };
    }

    public static AntSimulationData UpdatePheromoneAges(AntSimulationData simData, float deltaTime) {
        if(simData.Pheromones.Count > 0) {
            for (int i = simData.Pheromones.Count - 1; i >= 0; i--) {
                PheromoneData pheromone = simData.Pheromones[i];
                pheromone.age += pheromone.owner.pheromoneDecayRate * deltaTime;
                if (pheromone.age >= pheromone.owner.pheromoneMaxAge) simData.Pheromones.RemoveAt(i); // Remove expired pheromones
                else simData.Pheromones[i] = pheromone;
            }
        }
        return simData;
    }

    public static AntSimulationData ClearPheromonesInRangeAndPositionSensors(AntSimulationData simData) {
        // Clear the pheromonesInRange list for all ants
        for (int i = 0; i < simData.Ants.Count; i++) {
            AntData ant = simData.Ants[i];
            ant.pheromoneIndiciesInRange.Clear();
            ant.sensorPositions = CalculateSensorPositions(ant.baseData); // Calculate sensor positions for the ant
            simData.Ants[i] = ant;
        }
        return simData;
    }

    public static AntSimulationData UpdatePheromonesInRange(AntSimulationData simData) {
        // Iterate over each pheromone and check which ants' sensors are in range
        for (int pheromone = 0; pheromone < simData.Pheromones.Count; pheromone++) {
            bool relevant = false;
            for (int i = 0; i < simData.Ants.Count; i++) {
                AntData ant = simData.Ants[i];
                // Check if the pheromone is within the range of any of the ant's sensors
                foreach (float2 sensorPosition in ant.sensorPositions) {
                    float distance = Vector2.Distance(sensorPosition, simData.Pheromones[pheromone].position);
                    if (distance <= ant.baseData.sensorRadius) {
                        ant.pheromoneIndiciesInRange.Add(pheromone);
                        relevant = true;
                    }
                }
                simData.Ants[i] = ant;
            }
            if(!relevant && simData.PheromoneThreshold > simData.Pheromones.Count) {
                simData.Pheromones.RemoveAt(pheromone);
                pheromone--;
            }
        }
        return simData;
    }

    public static AntSimulationData UpdateAntEnergy(AntSimulationData simData) {
        for (int i = 0; i < simData.Ants.Count; i++) {
            AntData antData = simData.Ants[i]; // Get the struct
            AntJobData ant = antData.baseData;
            // Update ant energy
            ant.energy -= (ant.energyConsumptionRate * simData.SimulationTime);
            if (ant.energy <= 0) {
                simData.Ants.RemoveAt(i);
                i--;
            }
            else {
                antData.baseData = ant;
                simData.Ants[i] = antData;
            }
        }
        simData = AntCountCheck(simData);
        return simData;
    }

    public static AntSimulationData AntCountCheck(AntSimulationData simData) {
        int antCount = simData.Ants.Count;
        int colonySize = simData.ColonySize;
        if(antCount <= colonySize * 0.25) {
            int diff = colonySize - antCount;
            if (diff > 0) for (int i = 0; i < diff; i++) {
                simData.Ants.Add(new AntData {
                    baseData = AntSimulationCore.CreateOffspring(simData, math.min(antCount, diff / 10)),
                    pheromoneIndiciesInRange = new HashSet<int>(),
                    foodCollected = 0
                });
                antCount = simData.Ants.Count;
                diff = colonySize - antCount;
            }
        }
        return simData;
    }

    public static AntSimulationData UpdatePheromoneTimers(AntSimulationData simData) {
        for (int aIndex = 0; aIndex < simData.Ants.Count; aIndex++) {
            AntData antData = simData.Ants[aIndex];
            AntJobData ant = antData.baseData;
            ant.pheromoneTimer += simData.SimulationTime;
            antData.baseData = ant;
            simData.Ants[aIndex] = antData;
        }
        List<int> antsToUpdate = simData.Ants.Select((ant, index) => new { Ant = ant, Index = index }) // Pair each ant with its index
                                            .Where(pair => pair.Ant.baseData.pheromoneTimer >= pair.Ant.baseData.pheromoneIntervalValue)
                                            .Select(pair => pair.Index) // Select only the index
                                            .ToList();
        if(antsToUpdate.Count > 0) {
            simData = MergePheromones(simData, antsToUpdate);
            for (int i = 0; i < antsToUpdate.Count; i++) {
                int antIndex = antsToUpdate[i];
                AntData antData = simData.Ants[antIndex];
                AntJobData ant = antData.baseData;
                simData.Pheromones.Add(DepositPheromone(simData, ant));
                ant.pheromoneTimer = 0f;
                antData.baseData = ant;
                simData.Ants[antIndex] = antData;
            }
        }
        return simData;
    }

    public static AntSimulationData LocationCheck(AntSimulationData simData) {
        for (int i = 0; i < simData.Ants.Count; i++) {
            AntData antData = simData.Ants[i]; // Get the struct
            AntJobData ant = antData.baseData;
            // Check if ant found food
            if (ant.isCarryingFood == 0) {
                for (int j = 0; j < simData.FoodSources.Count; j++) {
                    FoodSource food = simData.FoodSources[j];
                    if (Vector2.Distance(ant.position, food.position) <= (float)food.size) {
                        ant.isCarryingFood = 1;
                        ant.rotation += 180f;
                        ant.rotation = NormaliseAngle(ant.rotation);
                        antData.foodCollected++;
                        // Update food source
                        FoodSource updatedFood = food;
                        updatedFood.mass -= 1f;
                        updatedFood.size = Mathf.Max((float)updatedFood.mass / 1000, 0.1f);
                        simData.FoodSources[j] = updatedFood;
                        break;
                    }
                }
            }
            // Check if ant reached nest
            else if (Vector2.Distance(ant.position, simData.NestPosition) <= simData.NestSize * 2) {
                // Ant delivered food!
                ant.isCarryingFood = 0;
                ant.energy = ant.maxEnergy; // Reward energy
                
                // Turn around to leave nest
                ant.rotation += 180f;
                ant.rotation = NormaliseAngle(ant.rotation);
            }
            antData.baseData = ant;
            simData.Ants[i] = antData;
        }
        return simData;
    }

    public static AntSimulationData MoveAnts(AntSimulationData simData) {
        for (int i = 0; i < simData.Ants.Count; i++) {
            AntData antData = simData.Ants[i]; // Get the struct
            AntJobData ant = antData.baseData;
            if (ant.position.x < -simData.WorldSize/2 || ant.position.x > simData.WorldSize/2 || ant.position.y < -simData.WorldSize/2 || ant.position.y > simData.WorldSize/2) {
                ant.rotation += 180 + random.NextFloat(float.Epsilon, 90);
                ant.rotation = NormaliseAngle(ant.rotation);
            }
            // Move ant
            ant.position += new float2(math.cos(ant.rotation * Mathf.Deg2Rad), math.sin(ant.rotation * Mathf.Deg2Rad)) * ant.speed * simData.SimulationTime;
            antData.baseData = ant;
            simData.Ants[i] = antData;
        }
        return simData;
    }

    public static AntSimulationData OrganiseLists(AntSimulationData simData) {
        simData.Ants = simData.Ants.OrderByDescending(ant => ant.foodCollected)
                                   .ThenByDescending(ant => ant.baseData.energy / ant.baseData.energyConsumptionRate)
                                   .ThenByDescending(ant => ant.baseData.speed)
                                   .ToList();
        simData.Pheromones = simData.Pheromones.OrderBy(p => p.owner.isCarryingFood)
                                    // .ThenBy(p => Vector2.Distance(p.position, simData.NestPosition)) // Sort by distance to nest
                                    .ThenByDescending(p => (p.age / p.owner.pheromoneMaxAge))
                                    .ToList();
        simData.FoodSources = simData.FoodSources.OrderByDescending(f => f.mass).ToList();
        return simData;
    }

    public static AntSimulationData Update(AntSimulationData simData, float deltaTime) {
        simData.SimulationTime = deltaTime;
        if (simData.Ants.Count > 0) {
            simData = UpdateAntEnergy(simData);
            simData = UpdatePheromoneTimers(simData);
            simData = LocationCheck(simData);
            simData = UpdateAntDirections(simData);
            if (simData.Pheromones.Count > 0) {
                simData = ClearPheromonesInRangeAndPositionSensors(simData);
                simData = UpdatePheromonesInRange(simData);
            }
            simData = UpdateAntRotations(simData);
            simData = MoveAnts(simData); // Move ants
        }
        simData = UpdateFoodSources(simData);
        return simData;
    }

    public static AntSimulationData UpdateFoodSources(AntSimulationData simData) {
        if (simData.FoodSources.Count > 0) {
            for (int i = simData.FoodSources.Count - 1; i >= 0; i--) if (simData.FoodSources[i].mass <= 0) simData.FoodSources.RemoveAt(i);
        }
        else simData.FoodSources.Add(CreateNewFoodSource(simData));
        return simData;
    }

    public static float2[] CalculateSensorPositions(AntJobData ant) {
        float2[] sensorPositions = new float2[3];
        float antDirection = ant.rotation * Mathf.Deg2Rad;
        // Front sensor
        sensorPositions[0] = ant.position + new float2(
            math.cos(antDirection) * ant.speed,
            math.sin(antDirection) * ant.speed
        );
        // Left sensor
        sensorPositions[1] = ant.position + new float2(
            math.cos(antDirection - ant.sensorAngle * Mathf.Deg2Rad) * ant.speed,
            math.sin(antDirection - ant.sensorAngle * Mathf.Deg2Rad) * ant.speed
        );
        // Right sensor
        sensorPositions[2] = ant.position + new float2(
            math.cos(antDirection + ant.sensorAngle * Mathf.Deg2Rad) * ant.speed,
            math.sin(antDirection + ant.sensorAngle * Mathf.Deg2Rad) * ant.speed
        );
        return sensorPositions;
    }
    
    public static AntSimulationData UpdateAntDirections(AntSimulationData simData) {
        for (int i = 0; i < simData.Ants.Count; i++) {
            AntData antData = simData.Ants[i]; // Get the struct
            AntJobData ant = antData.baseData;
            // Special case: If carrying food and near nest, steer toward nest
            if (ant.isCarryingFood == 1) {
                float distToNest = Vector2.Distance(ant.position, simData.NestPosition);
                if (distToNest <= simData.NestSize + ant.smellRadius) {
                    float2 dirToNest = simData.NestPosition - ant.position;
                    float targetAngle = Mathf.Atan2(dirToNest.y, dirToNest.x) * Mathf.Rad2Deg;
                    float angleDiff = Mathf.DeltaAngle(ant.rotation, targetAngle);
                    ant.rotation += Mathf.Clamp(angleDiff, -ant.rotationSpeed * simData.SimulationTime, ant.rotationSpeed * simData.SimulationTime);
                    ant.rotation = NormaliseAngle(ant.rotation);
                }
            }
            else {
            // Special case: If not carrying food and near food, steer toward food
                float closestFoodDistance = float.MaxValue;
                FoodSource? closestSource = null;
                for (int j = 0; j < simData.FoodSources.Count; j++) {
                    float distToFood = Vector2.Distance(ant.position, simData.FoodSources[j].position) - simData.FoodSources[j].size;
                    if (distToFood < closestFoodDistance) {
                        closestSource = simData.FoodSources[j];
                        closestFoodDistance = distToFood;
                    }
                }
                if (closestSource != null) {
                    if (closestFoodDistance <= ant.smellRadius + (closestSource.Value.size / 2)) {
                        Vector2 dirToFood = (float2)closestSource.Value.position - ant.position;
                        float targetAngle = Mathf.Atan2(dirToFood.y, dirToFood.x) * Mathf.Rad2Deg;
                        float angleDiff = Mathf.DeltaAngle(ant.rotation, targetAngle);
                        ant.rotation += Mathf.Clamp(angleDiff, -ant.rotationSpeed, ant.rotationSpeed);
                        ant.rotation = NormaliseAngle(ant.rotation);
                    }
                }
            }
            antData.baseData = ant;
            simData.Ants[i] = antData;
        }
        return simData;
    }

    public static AntSimulationData UpdateAntRotations(AntSimulationData simData) {
        int pheromoneCount = simData.Pheromones.Count;
        if (pheromoneCount > 0) {
            for (int i = 0; i < simData.Ants.Count; i++) {
                AntData antData = simData.Ants[i]; // Get the struct
                antData.sensorStrengths = new float[3] {0, 0, 0}; // Initialize sensor strengths to zero
                simData.Ants[i] = antData;
            }
            for (int pheromone = 0; pheromone < pheromoneCount; pheromone++) {
                for (int i = 0; i < simData.Ants.Count; i++) {
                    AntData antData = simData.Ants[i]; // Get the struct
                    if(antData.pheromoneIndiciesInRange.Contains(pheromone)) {
                        for (int j = 0; j < 3; j++) {
                            antData.sensorStrengths[j] += CalculatePheromoneConcentration(simData, antData, simData.Pheromones[pheromone], j);
                        }
                        for (int k = 0; k < 3; k++) {
                            antData.sensorStrengths[k] /= antData.pheromoneIndiciesInRange.Count;
                        }
                    }
                }
            }
            for (int i = 0; i < simData.Ants.Count; i++) {
                AntData antData = simData.Ants[i]; // Get the struct
                AntJobData ant = antData.baseData;
                // Normalize sensor strengths
                for (int j = 0; j < 3; j++) {
                    if (antData.sensorStrengths[j] != 0) antData.sensorStrengths[j] /= antData.pheromoneIndiciesInRange.Count;
                }
                bool leftSensorStrongest = (antData.sensorStrengths[1] >= antData.sensorStrengths[0] && antData.sensorStrengths[1] > antData.sensorStrengths[2]);
                bool rightSensorStrongest = (antData.sensorStrengths[2] >= antData.sensorStrengths[0] && antData.sensorStrengths[2] > antData.sensorStrengths[1]);
                // Calculate the direction of rotation based on sensor strengths
                float rotationDirection = (leftSensorStrongest) ? -1 : (rightSensorStrongest) ? 1 : 0;
                float rotationChange = (rotationDirection * ant.rotationSpeed * simData.SimulationTime);
                // Apply the rotation
                ant.rotation += rotationChange;
                ant.rotation = NormaliseAngle(ant.rotation);
                antData.baseData = ant;
                simData.Ants[i] = antData;
            }
        }
        return simData;
    }

    private static float CalculatePheromoneConcentration (AntSimulationData simData, AntData antData, PheromoneData pheromone, int SensorIndex) {
        AntJobData ant = antData.baseData;
        float sensorDistanceToPheromone = Vector2.Distance(pheromone.position, antData.sensorPositions[SensorIndex]);
        float totalStrength = 0;
        if (sensorDistanceToPheromone <= ant.sensorRadius) {
            int antHasFood = ant.isCarryingFood;
            int pheromoneOwnerHadFood = pheromone.owner.isCarryingFood;
            float distanceHome = Vector2.Distance(pheromone.position, simData.NestPosition);
            float antDistanceHome = Vector2.Distance(ant.position, simData.NestPosition);
            float distanceRatio = (antHasFood == 1) ? antDistanceHome / (distanceHome + float.Epsilon) : distanceHome / (antDistanceHome + float.Epsilon);

            float directionStrength = (antHasFood == 1) ? ((distanceHome < antDistanceHome) ? distanceRatio * (simData.movementSettings.directionSensitivity * 100) : -distanceRatio) :
                                                          ((distanceHome > antDistanceHome) ? distanceRatio * (simData.movementSettings.directionSensitivity * 100) : -distanceRatio);
            float homeStrength = (antDistanceHome / (distanceHome + float.Epsilon)) * (simData.movementSettings.distanceHomeSensitivity * 100);
            float pheromoneMaturity = ((pheromone.age + float.Epsilon) / (pheromone.owner.pheromoneMaxAge + float.Epsilon)) * (simData.movementSettings.pheromoneAgeSensitivity * 100);
            float typestrength = (antHasFood == 1) ? ((pheromoneOwnerHadFood == 0) ? (simData.movementSettings.pheromoneTypeSensitivity * 100) : 1) :
                                                     ((pheromoneOwnerHadFood == 1) ? (simData.movementSettings.pheromoneTypeSensitivity * 100) : 1);
            float sensorDistanceStrength = (sensorDistanceToPheromone / ant.sensorRadius) * (simData.movementSettings.sensorDistanceSensitivity * 100);

            float distanceStrength = (directionStrength) * (pheromoneMaturity / 1);
            totalStrength = ((sensorDistanceStrength) * (distanceStrength)) * typestrength;
        }
        return totalStrength;
    }

    public static float NormaliseAngle(float Angle) {
        float newAngle = Angle % 360f; // Normalize rotation
        if (newAngle < 0) newAngle += 360f; // Between 0 and 359
        return newAngle;
    }

    public static AntSimulationData MergePheromones(AntSimulationData simData, List<int> antIndices) {
        int pheromoneCount = simData.Pheromones.Count;
        if (pheromoneCount > 0) {
            for(int i = pheromoneCount - 1; i >= 0; i--) {
                PheromoneData pheromone = simData.Pheromones[i];
                int antCount = antIndices.Count;
                if(antCount > 0) {
                    for (int a = 0; a < antCount; a++) {
                        AntJobData ant = simData.Ants[antIndices[a]].baseData;
                        int pointsToFood = ant.isCarryingFood; // Points to food if ant is carrying food
                        if (Vector2.Distance(pheromone.position, ant.position) <= ant.pheromoneMergeRadius) {
                            if(pheromone.owner.isCarryingFood == pointsToFood || pointsToFood == 1) {
                                simData.Pheromones.RemoveAt(i);
                                i--; // Adjust index after removal
                                break;
                            }
                        }
                    }
                }
            }
        }
        return simData;
    }

    public static PheromoneData DepositPheromone(AntSimulationData simData, AntJobData ant) {
        return new PheromoneData {
            position = ant.position,
            owner = ant,
            age = 0f
        };
    }
    
    public static FoodSource CreateNewFoodSource(AntSimulationData simData) {
        Vector2 pos;
        do {
            pos = new Vector2(
                UnityEngine.Random.Range(-simData.WorldSize/2, simData.WorldSize/2),
                UnityEngine.Random.Range(-simData.WorldSize/2, simData.WorldSize/2)
            );
        } while (Vector2.Distance(pos, simData.NestPosition) < simData.WorldSize * 0.05f);
        float massa = UnityEngine.Random.Range(500f, 20000f);
        return new FoodSource {
            position = pos,
            mass = massa,
            size = Mathf.Max(massa / 1000f, 0.01f)
        };
    }
}

public interface ISimulationRenderer {
    void Initialize();
    void Draw(AntSimulationData data, bool drawPheromones);
    void CleanUp();
}
#region Structs
[System.Serializable]
public struct AntSimulationData {
    public float SimulationTime;
    public int ColonySize;
    public float2 NestPosition;
    public float NestSize;
    public float WorldSize;
    public int PheromoneThreshold;
    public List<AntData> Ants;
    public List<PheromoneData> Pheromones;
    public List<FoodSource> FoodSources;
    public MovementParamaters movementSettings;
}

[System.Serializable]
public struct AntJobData {
    public float2 position;
    public float rotation;
    public int isCarryingFood;

    public float pheromoneTimer;
    public float pheromoneIntervalValue;
    public float pheromoneDecayRate;
    public float pheromoneMaxAge;
    public float pheromoneMergeRadius;

    public float sensorRadius;
    public float sensorAngle;

    public float smellRadius;

    public float speed;
    public float rotationSpeed;

    public float energy;
    public float energyConsumptionRate;
    public float maxEnergy;

}

[System.Serializable]
public struct AntData {
    public AntJobData baseData;
    public HashSet<int> pheromoneIndiciesInRange;
    public float2[] sensorPositions;
    public float[] sensorStrengths;
    public int foodCollected;
}

[System.Serializable]
public struct PheromoneData {
    public float2 position;
    public AntJobData owner;
    public float age;
}

[System.Serializable]
public struct FoodSource {
    public Vector2 position;
    public float mass;
    public float size;
}

[System.Serializable]
public struct MovementParamaters {
    public float pheromoneAgeSensitivity;
    public float pheromoneTypeSensitivity;
    public float distanceHomeSensitivity;
    public float directionSensitivity;
    public float sensorDistanceSensitivity;
}

[BurstCompile]
public struct CreateAntsJob : IJobParallelFor {
    [ReadOnly] public float2 NestPosition;
    [ReadOnly] public float NestSize;
    [WriteOnly] public NativeArray<AntJobData> Ants;

    public void Execute(int index) {
        Ants[index] = AntSimulationCore.CreateAnt(NestPosition, NestSize);
    }
}

[BurstCompile]
public struct UpdatePheromoneAgeJob : IJobParallelFor {
    public NativeArray<PheromoneData> PheromoneData;
    [ReadOnly] public float DeltaTime;
    
    public NativeQueue<int>.ParallelWriter RemoveQueueWriter;

    public void Execute(int index) {
        var pheromone = PheromoneData[index];
        pheromone.age += pheromone.owner.pheromoneDecayRate * DeltaTime;
        if (pheromone.age >= pheromone.owner.pheromoneMaxAge) RemoveQueueWriter.Enqueue(index);
        PheromoneData[index] = pheromone;
    }
}

[BurstCompile]
public struct ProcessPheromonesJob : IJob {
    [ReadOnly] public NativeArray<AntJobData> Ants;
    [ReadOnly] public NativeArray<PheromoneData> Pheromones;
    [ReadOnly] public NativeArray<float2> SensorPositions; // Flattened array: 3 sensors per ant
    public UnsafeParallelMultiHashMap<int, int> PheromonesInRangeMap; // Each job has its own map
    public int StartIndex;
    public int EndIndex;

    public void Execute() {
        for (int pheromoneIndex = StartIndex; pheromoneIndex < EndIndex; pheromoneIndex++) {
            bool relevant = false;
            PheromoneData pheromone = Pheromones[pheromoneIndex];
            for (int antIndex = 0; antIndex < Ants.Length; antIndex++) {
                AntJobData ant = Ants[antIndex];
                // Check each of the ant's 3 sensors
                for (int s = 0; s < 3; s++) {
                    int sensorIdx = (antIndex * 3) + s;
                    float2 sensorPos = SensorPositions[sensorIdx];
                    float dist = math.distance(sensorPos, pheromone.position);
                    if (dist <= ant.sensorRadius) {
                        PheromonesInRangeMap.Add(antIndex, pheromoneIndex);
                        relevant = true;
                        break; // No need to check other sensors
                    }
                }
            }
            if (!relevant) {
                // Mark this pheromone for removal
                PheromonesInRangeMap.Add(-1, pheromoneIndex);
            }
        }
    }
}
#endregion
