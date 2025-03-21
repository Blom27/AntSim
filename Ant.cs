using UnityEngine;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;
using Unity.Burst;
using System.Collections.Generic;
using Unity.Collections.LowLevel.Unsafe;
using System.Linq;

public class AntSimulationController : MonoBehaviour {
    public enum RenderMode { Gizmos/*, Sprites, ThreeD*/ }
    
    [Header("Simulation Settings")]
    public int colonySize = 50;
    public float2 nestPosition = Vector2.zero;
    public float worldSize = 50f;
    public RenderMode currentRenderMode = RenderMode.Gizmos;
    public ISimulationRenderer _activeRenderer;
    
    [Header("Data")]
    public AntSimulationData _simulationData;

    private NativeArray<PheromoneData> _pheromoneJobData;

    void Start() {
        _simulationData = AntSimulationCore.Initialize(_simulationData, colonySize, nestPosition, 3, worldSize, 0.5f);
        if(colonySize > 256) CreateAntsWithJobSystem(colonySize);
        else _simulationData = AntSimulationCore.CreateAnts(_simulationData);
        SwitchRenderer(currentRenderMode);
    }

    void Update() {
        if (_simulationData.Ants.Count > 0) {
            if (Input.GetMouseButtonDown(0)) _simulationData.FoodSources.Add(AntSimulationCore.CreateNewFoodSource(_simulationData));
            if(_simulationData.Pheromones.Count > 64){
                _simulationData = AntSimulationCore.UpdateAntEnergy(_simulationData);
                UpdatePheromonesWithJobSystem(Time.deltaTime);
                _simulationData = AntSimulationCore.UpdatePheromoneTimers(_simulationData);
                _simulationData = AntSimulationCore.LocationCheck(_simulationData);
                _simulationData = AntSimulationCore.ClearPheromonesInRangeAndPositionSensors(_simulationData);
                UpdatePheromonesInRangeWithJobSystem();
                _simulationData = AntSimulationCore.UpdateAntRotations(_simulationData);
                _simulationData = AntSimulationCore.UpdateAntDirections(_simulationData);
                _simulationData = AntSimulationCore.MoveAnts(_simulationData); // Move ants
                _simulationData = AntSimulationCore.UpdateFoodSources(_simulationData);
            }
            else _simulationData = AntSimulationCore.Update(_simulationData, Time.deltaTime);
        }
        if(currentRenderMode != RenderMode.Gizmos) _activeRenderer.Draw(_simulationData);
    }

    void OnDrawGizmos() {
        if(currentRenderMode == RenderMode.Gizmos && Application.isPlaying) {
            _activeRenderer.Draw(_simulationData);
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
            Ants = antJobData,
            Seed = (uint)System.DateTime.Now.Millisecond
        };
        
        // Execute the job
        JobHandle jobHandle = createAntsJob.Schedule(colonySize, 64);
        jobHandle.Complete();
        
        // Transfer results to AntSimulationCore
        _simulationData.Ants = antJobData.Select(ajd => new AntData {
            baseData = ajd,
            pheromonesInRange = new List<PheromoneData>()
        }).ToList();

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
        UnsafeParallelMultiHashMap<int, PheromoneData> pheromonesInRangeMap = new UnsafeParallelMultiHashMap<int, PheromoneData>(pheromoneCount * 3, Allocator.TempJob);

        NativeParallelHashSet<int> antsProcessedMap = new NativeParallelHashSet<int>(antCount + 1, Allocator.TempJob);

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

        // Set up and schedule the job
        var updatePheromonesInRangeJob = new UpdatePheromonesInRangeJob {
            Ants = antJobData,
            Pheromones = pheromoneJobData,
            SensorPositions = sensorPositions,
            PheromonesInRangeMap = pheromonesInRangeMap.AsParallelWriter(),
            AntsProcessedMap = antsProcessedMap.AsParallelWriter()
        };

        JobHandle jobHandle = updatePheromonesInRangeJob.Schedule(pheromoneCount, 64);
        jobHandle.Complete();

        // Transfer results back to the simulation data
        for (int i = 0; i < antCount; i++) {
            _simulationData.Ants[i].pheromonesInRange.Clear();
        }

        foreach (var pair in pheromonesInRangeMap) {
            if (!_simulationData.Ants[pair.Key].pheromonesInRange.Contains(pair.Value)) {
                _simulationData.Ants[pair.Key].pheromonesInRange.Add(pair.Value);
            }
        }

        // Dispose native containers
        antJobData.Dispose();
        pheromoneJobData.Dispose();
        sensorPositions.Dispose();
        pheromonesInRangeMap.Dispose();
        antsProcessedMap.Dispose();
    }
}

[System.Serializable]
public class GizmoRenderer : ISimulationRenderer {
    private Color _antColor = Color.white;
    private Color _foodColor = Color.green;

    public void Initialize() { /* Not needed for Gizmos */ }

    public void Draw(AntSimulationData data) {
        DrawPheromones(data.Pheromones);
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
    public static AntSimulationData currentState;

    public static Unity.Mathematics.Random random = new Unity.Mathematics.Random((uint)System.DateTime.Now.Millisecond);

    public static AntSimulationData Initialize(AntSimulationData simData, int colonySize, Vector2 nestPos, int foodSourceCount, float WorldSize, float NestSize) {
        simData = new AntSimulationData {
            SimulationTime = 0f,
            ColonySize = colonySize,
            NestPosition = nestPos,
            NestSize = NestSize,
            WorldSize = WorldSize,
            Ants = new List<AntData>(),
            Pheromones = new List<PheromoneData>(),
            FoodSources = new List<FoodSource>()
        };
        // Create food sources
        for(int i = 0; i < foodSourceCount; i++) simData.FoodSources.Add(CreateNewFoodSource(simData));
        return simData;
    }

    public static AntSimulationData CreateAnts(AntSimulationData simData) {
        for (int i = 0; i < simData.ColonySize; i++) {
            AntJobData ant = CreateAnt(simData.NestPosition, simData.NestSize);
            simData.Ants.Add(new AntData {
                baseData = ant,
                pheromonesInRange = new List<PheromoneData>(),
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

        float PheromoneIntervalValue = random.NextFloat(0f, 1f);
        float PheromoneDecayRate = random.NextFloat(0f, 0.1f);
        float PheromoneMaxAge = PheromoneDecayRate / PheromoneIntervalValue;
        float PheromoneMergeRadius = random.NextFloat(0f, 1f) / 10;
        float Speed = random.NextFloat(0f, 1f) * 2;
        float RotationSpeed = random.NextFloat(0f, 360f) / 2;
        float MaxEnergy = ((PheromoneDecayRate + PheromoneMergeRadius + Speed + (RotationSpeed / 360)) / (PheromoneIntervalValue + PheromoneMaxAge)) * 15;
        float SensorAngle = random.NextFloat(30f, 60f); // Randomize sensor angle

        return new AntJobData {
            position = position,
            rotation = rotation,
            isCarryingFood = 0,

            pheromoneTimer = 0f,
            pheromoneIntervalValue = PheromoneIntervalValue,
            pheromoneDecayRate = PheromoneDecayRate,
            pheromoneMaxAge = PheromoneMaxAge,
            pheromoneMergeRadius = PheromoneMergeRadius,

            sensorRadius = random.NextFloat(0.5f, 1f) * 2,
            sensorAngle = SensorAngle,
            smellRadius = random.NextFloat(1f, 2f) * 2,

            speed = Speed,
            rotationSpeed = RotationSpeed,
            maxEnergy = MaxEnergy,
            energy = MaxEnergy,
            energyConsumptionRate = (PheromoneMaxAge * Speed * (RotationSpeed / 180)) / 100

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
            ant.pheromonesInRange.Clear();
            ant.sensorPositions = CalculateSensorPositions(ant.baseData); // Calculate sensor positions for the ant
            simData.Ants[i] = ant;
        }
        return simData;
    }

    public static AntSimulationData UpdatePheromonesInRange(AntSimulationData simData) {
        // Iterate over each pheromone and check which ants' sensors are in range
        for (int pheromone = 0; pheromone < simData.Pheromones.Count; pheromone++) {
            for (int i = 0; i < simData.Ants.Count; i++) {
                AntData ant = simData.Ants[i];
                // Check if the pheromone is within the range of any of the ant's sensors
                foreach (float2 sensorPosition in ant.sensorPositions) {
                    float distance = Vector2.Distance(sensorPosition, simData.Pheromones[pheromone].position);
                    if (distance <= ant.baseData.sensorRadius) ant.pheromonesInRange.Add(simData.Pheromones[pheromone]);
                }
                simData.Ants[i] = ant;
            }
        }
        return simData;
    }

    public static AntSimulationData UpdateAntEnergy(AntSimulationData simData) {
        for (int i = 0; i < simData.Ants.Count; i++) {
            AntData antData = simData.Ants[i]; // Get the struct
            AntJobData ant = antData.baseData;
            // Update ant energy
            ant.energy -= ant.energyConsumptionRate * simData.SimulationTime;
            if (ant.energy <= 0) {
                simData.Ants.RemoveAt(i);
                i--;
            }
            else {
                antData.baseData = ant;
                simData.Ants[i] = antData;
            }
        }
        return simData;
    }

    public static AntSimulationData UpdatePheromoneTimers(AntSimulationData simData) {
        for (int i = 0; i < simData.Ants.Count; i++) {
            AntData antData = simData.Ants[i]; // Get the struct
            AntJobData ant = antData.baseData;
            // Update pheromone timer and deposit if needed
            ant.pheromoneTimer += simData.SimulationTime;
            if (ant.pheromoneTimer >= ant.pheromoneIntervalValue) {
                ant = DepositPheromone(simData, ant);
                ant.pheromoneTimer = 0f;
            }
            antData.baseData = ant;
            simData.Ants[i] = antData;
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
            else if (Vector2.Distance(ant.position, simData.NestPosition) <= simData.NestSize) {
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
            // Move ant
            ant.position += new float2(math.cos(ant.rotation * Mathf.Deg2Rad), math.sin(ant.rotation * Mathf.Deg2Rad)) * ant.speed * simData.SimulationTime;
            if (ant.position.x < -simData.WorldSize/2 || ant.position.x > simData.WorldSize/2 || ant.position.y < -simData.WorldSize/2 || ant.position.y > simData.WorldSize/2) {
                ant.position = simData.NestPosition;
                ant.rotation += 180f;
                ant.rotation = NormaliseAngle(ant.rotation);
            }
            antData.baseData = ant;
            simData.Ants[i] = antData;
        }
        return simData;
    }

    public static AntSimulationData Update(AntSimulationData simData, float deltaTime) {
        simData.SimulationTime = deltaTime;
        if (simData.Ants.Count > 0) {
            simData = UpdateAntEnergy(simData);
            simData = UpdatePheromoneTimers(simData);
            simData = LocationCheck(simData);
            if (simData.Pheromones.Count > 0) {
                simData = ClearPheromonesInRangeAndPositionSensors(simData);
                simData = UpdatePheromonesInRange(simData);
            }
            simData = UpdateAntRotations(simData);
            simData = UpdateAntDirections(simData);
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
            math.cos(antDirection) * ant.sensorRadius,
            math.sin(antDirection) * ant.sensorRadius
        );
        // Left sensor
        sensorPositions[1] = ant.position + new float2(
            math.cos(antDirection - ant.sensorAngle * Mathf.Deg2Rad) * ant.sensorRadius,
            math.sin(antDirection - ant.sensorAngle * Mathf.Deg2Rad) * ant.sensorRadius
        );
        // Right sensor
        sensorPositions[2] = ant.position + new float2(
            math.cos(antDirection + ant.sensorAngle * Mathf.Deg2Rad) * ant.sensorRadius,
            math.sin(antDirection + ant.sensorAngle * Mathf.Deg2Rad) * ant.sensorRadius
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
                    if (closestFoodDistance <= ant.smellRadius + closestSource.Value.size) {
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
        for (int i = 0; i < simData.Ants.Count; i++) {
            AntData antData = simData.Ants[i]; // Get the struct
            AntJobData ant = antData.baseData;
            antData.sensorStrengths = new float[3] {0, 0, 0}; // Initialize sensor strengths to zero
            if (antData.pheromonesInRange.Count > 0) {
                foreach (PheromoneData pheromone in antData.pheromonesInRange) {
                    for (int j = 0; j < antData.sensorPositions.Length; j++) {
                        antData.sensorStrengths[j] += CalculatePheromoneConcentration(simData, antData, pheromone, j);
                    }
                }
                // Normalize sensor strengths
                for (int j = 0; j < 3; j++) {
                    if (antData.sensorStrengths[j] != 0) antData.sensorStrengths[j] /= antData.pheromonesInRange.Count;
                }
                bool leftSensorStrongest = (antData.sensorStrengths[1] >= antData.sensorStrengths[0] && antData.sensorStrengths[1] > antData.sensorStrengths[2]);
                bool rightSensorStrongest = (antData.sensorStrengths[2] >= antData.sensorStrengths[0] && antData.sensorStrengths[2] > antData.sensorStrengths[1]);
                // Calculate the direction of rotation based on sensor strengths
                float rotationDirection = (leftSensorStrongest) ? -1 : (rightSensorStrongest) ? 1 : 0;
                float rotationChange = (rotationDirection * ant.rotationSpeed * simData.SimulationTime);
                // Apply the rotation
                ant.rotation += rotationChange;
                ant.rotation = NormaliseAngle(ant.rotation);
            }
            antData.baseData = ant;
            simData.Ants[i] = antData;
        }
        return simData;
    }

    private static float CalculatePheromoneConcentration (AntSimulationData simData, AntData antData, PheromoneData pheromone, int SensorIndex) {
        AntJobData ant = antData.baseData;
        float sensorDistanceToPheromone = Vector2.Distance(pheromone.position, antData.sensorPositions[SensorIndex]);
        bool relevantPheromone = (ant.isCarryingFood == 1) ? true : (pheromone.owner.isCarryingFood == 1);
        float totalStrength = 0;
        if (sensorDistanceToPheromone <= ant.sensorRadius && relevantPheromone) {
            float distanceHome = Vector2.Distance(pheromone.position, simData.NestPosition);
            float antDistanceHome = Vector2.Distance(ant.position, simData.NestPosition);
            float distanceRatio = (ant.isCarryingFood == 1) ? antDistanceHome / (distanceHome + float.Epsilon) : distanceHome / (antDistanceHome + float.Epsilon);
            float directionStrength = (ant.isCarryingFood == 1) ? ((distanceHome < antDistanceHome) ? distanceRatio : -distanceRatio) :
                                                                ((distanceHome > antDistanceHome) ? distanceRatio : -distanceRatio);
            float distanceStrength = distanceHome / directionStrength;
            float sensorDistanceStrength = sensorDistanceToPheromone / ant.sensorRadius;
            totalStrength = (distanceStrength * sensorDistanceStrength);
        }
        return totalStrength;
    }

    public static float NormaliseAngle(float Angle) {
        float newAngle = Angle % 360f; // Normalize rotation
        if (newAngle < 0) newAngle += 360f; // Between 0 and 359
        return newAngle;
    }

    public static AntJobData DepositPheromone(AntSimulationData simData, AntJobData ant) {
        int pointsToFood = ant.isCarryingFood; // Points to food if ant is carrying food
        for(int i = 0; i < simData.Pheromones.Count; i++) {
            PheromoneData pheromone = simData.Pheromones[i];
            if (Vector2.Distance(pheromone.position, ant.position) <= ant.pheromoneMergeRadius) {
                if(pointsToFood == pheromone.owner.isCarryingFood) {
                    simData.Pheromones.RemoveAt(i);
                    i--; // Adjust index after removal
                }
            }
        }
        simData.Pheromones.Add(new PheromoneData {
            position = ant.position,
            owner = ant,
            age = 0f
        });
        return ant;
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
    void Draw(AntSimulationData data);
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
    public List<AntData> Ants;
    public List<PheromoneData> Pheromones;
    public List<FoodSource> FoodSources;
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
    public List<PheromoneData> pheromonesInRange;
    public float2[] sensorPositions;
    public float[] sensorStrengths;
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

[BurstCompile]
public struct CreateAntsJob : IJobParallelFor {
    [ReadOnly] public float2 NestPosition;
    [ReadOnly] public float NestSize;
    public uint Seed;
    [WriteOnly] public NativeArray<AntJobData> Ants;

    public void Execute(int index) {
        var random = new Unity.Mathematics.Random(Seed + (uint)index);
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
public struct UpdatePheromonesInRangeJob : IJobParallelFor {
    [ReadOnly] public NativeArray<AntJobData> Ants;
    [ReadOnly] public NativeArray<PheromoneData> Pheromones;
    [ReadOnly] public NativeArray<float2> SensorPositions; // Flattened array: 3 sensors per ant
    [WriteOnly] public UnsafeParallelMultiHashMap<int, PheromoneData>.ParallelWriter PheromonesInRangeMap;
    public NativeParallelHashSet<int>.ParallelWriter AntsProcessedMap;

    public void Execute(int index) {
        PheromoneData pheromone = Pheromones[index];
        for (int antIndex = 0; antIndex < Ants.Length; antIndex++) {
            AntJobData ant = Ants[antIndex];
            // Check each of the ant's 3 sensors
            for (int s = 0; s < 3; s++) {
                int sensorIdx = (antIndex * 3) + s;
                float2 sensorPos = SensorPositions[sensorIdx];
                float dist = math.distance(sensorPos, pheromone.position);
                if (dist <= ant.sensorRadius) {
                    
                    // Check if the pheromone has already been added for this ant
                    if(AntsProcessedMap.Add(antIndex)) {
                        PheromonesInRangeMap.Add(antIndex, pheromone);
                        break; // No need to check other sensors
                    }
                }
            }
        }
    }
}
#endregion
