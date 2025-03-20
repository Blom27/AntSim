using UnityEngine;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;
using Unity.Burst;
using System.Collections.Generic;
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
        AntSimulationCore.Initialize(colonySize, nestPosition, 3, worldSize, 0.5f);
        if(colonySize > 256) CreateAntsWithJobSystem(colonySize);
        else AntSimulationCore.CreateAnts(colonySize);
        SwitchRenderer(currentRenderMode);
    }

    void Update() {
        if (Input.GetMouseButtonDown(0)) AntSimulationCore.foodSources.Add(AntSimulationCore.CreateNewFoodSource());
        if(AntSimulationCore.pheromones.Count > 256) UpdatePheromonesWithJobSystem(Time.deltaTime);
        else AntSimulationCore.UpdatePheromoneAges(Time.deltaTime);
        AntSimulationCore.Update(Time.deltaTime);
        UpdateSimulationData();
        
        if(currentRenderMode != RenderMode.Gizmos) _activeRenderer.Draw(_simulationData);
    }

    void OnDrawGizmos() {
        if(currentRenderMode == RenderMode.Gizmos && Application.isPlaying) {
            _activeRenderer.Draw(_simulationData);
        }
    }

    void UpdateSimulationData() {
        _simulationData = new AntSimulationData {
            SimulationTime = AntSimulationCore.simulationTime,
            NestPosition = AntSimulationCore.nestPosition,
            NestSize = AntSimulationCore.nestSize,
            WorldSize = AntSimulationCore.worldSize,
            Ants = AntSimulationCore.ants,
            Pheromones = AntSimulationCore.pheromones,
            FoodSources = AntSimulationCore.foodSources
        };
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
            NestSize = (float)AntSimulationCore.nestSize,
            Ants = antJobData,
            Seed = (uint)System.DateTime.Now.Millisecond
        };
        
        // Execute the job
        JobHandle jobHandle = createAntsJob.Schedule(colonySize, 64);
        jobHandle.Complete();
        
        // Transfer results to AntSimulationCore
        AntSimulationCore.ants = antJobData.Select(ajd => new AntData {
            baseData = ajd,
            pheromonesInRange = new List<PheromoneData>()
        }).ToList();

        antJobData.Dispose();
    }

    private void UpdatePheromonesWithJobSystem(float deltaTime) {
        int pheromoneCount = AntSimulationCore.pheromones.Count;
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
                    position = new float2(AntSimulationCore.pheromones[i].position.x, 
                                        AntSimulationCore.pheromones[i].position.y),
                    owner = AntSimulationCore.pheromones[i].owner,
                    age = AntSimulationCore.pheromones[i].age
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
            AntSimulationCore.pheromones = newPheromones;
            indicesToRemove.Dispose();
        }
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
        
        // Add this required field
        public NativeQueue<int>.ParallelWriter RemoveQueueWriter;

        public void Execute(int index) {
            var pheromone = PheromoneData[index];
            pheromone.age += pheromone.owner.pheromoneDecayRate * DeltaTime;
            if (pheromone.age >= pheromone.owner.pheromoneMaxAge) RemoveQueueWriter.Enqueue(index);
            PheromoneData[index] = pheromone;
        }
    }
}

[System.Serializable]
public struct AntSimulationData {
    public float SimulationTime;
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

public static class AntSimulationCore {
    public static float simulationTime = 0f;
    public static float2 nestPosition;
    public static float nestSize;
    public static float worldSize = 50f;
    public static List<AntData> ants = new List<AntData>();
    public static List<PheromoneData> pheromones = new List<PheromoneData>();
    public static List<FoodSource> foodSources = new List<FoodSource>();
    public static Unity.Mathematics.Random random = new Unity.Mathematics.Random((uint)System.DateTime.Now.Millisecond);

    public static void Initialize(int colonySize, Vector2 nestPos, int foodSourceCount, float WorldSize, float NestSize) {
        nestPosition = nestPos;
        ants = new List<AntData>();
        pheromones = new List<PheromoneData>();
        foodSources = new List<FoodSource>();
        simulationTime = 0f;
        worldSize = WorldSize;
        nestSize = NestSize;
        
        // Create food sources
        for(int i = 0; i < foodSourceCount; i++) foodSources.Add(CreateNewFoodSource());
    }

    public static void CreateAnts(int ColonySize) {
        for (int i = 0; i < ColonySize; i++) {
            AntJobData ant = CreateAnt(nestPosition, nestSize);
            ants.Add(new AntData {
                baseData = ant,
                pheromonesInRange = new List<PheromoneData>(),
                sensorPositions = new float2[3],
                sensorStrengths = new float[3]
            });
        }
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

    public static void UpdatePheromoneAges(float deltaTime) {
        if(pheromones.Count > 0) {
            for (int i = pheromones.Count - 1; i >= 0; i--) {
                PheromoneData pheromone = pheromones[i];
                pheromone.age += pheromone.owner.pheromoneDecayRate * deltaTime;
                if (pheromone.age >= pheromone.owner.pheromoneMaxAge) pheromones.RemoveAt(i); // Remove expired pheromones
                else pheromones[i] = pheromone;
            }
        }
    }

    public static void UpdatePheromonesInRange() {
        // Clear the pheromonesInRange list for all ants
        foreach (AntData antData in ants) antData.pheromonesInRange.Clear();
        // Iterate over each pheromone and check which ants' sensors are in range
        for (int pheromone = 0; pheromone < pheromones.Count; pheromone++) {
            for (int i = 0; i < ants.Count; i++) {
                AntData ant = ants[i];
                // Calculate sensor positions for the ant
                float2[] sensorPositions = CalculateSensorPositions(ant.baseData);
                // Check if the pheromone is within the range of any of the ant's sensors
                foreach (float2 sensorPosition in sensorPositions) {
                    float distance = Vector2.Distance(sensorPosition, pheromones[pheromone].position);
                    if (distance <= ant.baseData.sensorRadius) ant.pheromonesInRange.Add(pheromones[pheromone]);
                }
                ant.sensorPositions = sensorPositions;
                ants[i] = ant;
            }
        }
    }

    public static void UpdateAnts(float deltaTime) {
        for (int i = 0; i < ants.Count; i++) {
            AntData antData = ants[i]; // Get the struct
            AntJobData ant = antData.baseData;
            // Update ant energy
            ant.energy -= ant.energyConsumptionRate * deltaTime;
            if (ant.energy <= 0) {
                ants.RemoveAt(i);
                i--;
                continue;
            }
            // Update pheromone timer and deposit if needed
            ant.pheromoneTimer += deltaTime;
            if (ant.pheromoneTimer >= ant.pheromoneIntervalValue) {
                DepositPheromone(ant);
                ant.pheromoneTimer = 0f;
            }
            // Check if ant found food
            if (ant.isCarryingFood == 0) {
                for (int j = 0; j < foodSources.Count; j++) {
                    FoodSource food = foodSources[j];
                    if (Vector2.Distance(ant.position, food.position) <= (float)food.size) {
                        ant.isCarryingFood = 1;
                        ant.rotation += 180f;
                        ant.rotation = NormaliseAngle(ant.rotation);
                        // Update food source
                        FoodSource updatedFood = food;
                        updatedFood.mass -= 1f;
                        updatedFood.size = Mathf.Max((float)updatedFood.mass / 1000, 0.1f);
                        foodSources[j] = updatedFood;
                        break;
                    }
                }
            }
            // Check if ant reached nest
            else if (Vector2.Distance(ant.position, nestPosition) <= nestSize) {
                // Ant delivered food!
                ant.isCarryingFood = 0;
                ant.energy = ant.maxEnergy; // Reward energy
                
                // Turn around to leave nest
                ant.rotation += 180f;
                ant.rotation = NormaliseAngle(ant.rotation);
            }
            // Update ant direction based on pheromones
            UpdateAntDirection(ref antData, deltaTime);
            if (antData.pheromonesInRange.Count > 0) antData = UpdateAntRotation(antData);
            // Move ant
            ant.position += new float2(math.cos(ant.rotation * Mathf.Deg2Rad), math.sin(ant.rotation * Mathf.Deg2Rad)) * ant.speed * deltaTime;
            if (ant.position.x < -worldSize/2 || ant.position.x > worldSize/2 || ant.position.y < -worldSize/2 || ant.position.y > worldSize/2) {
                ant.position = nestPosition;
                ant.rotation += 180f;
                ant.rotation = NormaliseAngle(ant.rotation);
            }
            antData.baseData = ant;
            ants[i] = antData;
        }
    }

    public static void Update(float deltaTime) {
        simulationTime += deltaTime;
        if (pheromones.Count > 0) UpdatePheromonesInRange();
        if (ants.Count > 0) UpdateAnts(deltaTime);
        UpdateFoodSources();
    }

    public static void UpdateFoodSources() {
        if (foodSources.Count > 0) {
            for (int i = foodSources.Count - 1; i >= 0; i--) if (foodSources[i].mass <= 0) foodSources.RemoveAt(i);
        }
        else foodSources.Add(CreateNewFoodSource());
    }

    private static float2[] CalculateSensorPositions(AntJobData ant) {
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
    
    private static void UpdateAntDirection(ref AntData antData, float deltaTime) {
        AntJobData ant = antData.baseData;
        // Special case: If carrying food and near nest, steer toward nest
        if (ant.isCarryingFood == 1) {
            float distToNest = Vector2.Distance(ant.position, nestPosition);
            if (distToNest < nestSize + ant.smellRadius) {
                float2 dirToNest = (float2)nestPosition - ant.position;
                float targetAngle = Mathf.Atan2(dirToNest.y, dirToNest.x) * Mathf.Rad2Deg;
                float angleDiff = Mathf.DeltaAngle(ant.rotation, targetAngle);
                ant.rotation += Mathf.Clamp(angleDiff, -ant.rotationSpeed, ant.rotationSpeed);
                ant.rotation = NormaliseAngle(ant.rotation);
                return;
            }
        }
        else {
        // Special case: If not carrying food and near food, steer toward food
            float closestFoodDistance = float.MaxValue;
            FoodSource? closestSource = null;
            for (int i = 0; i < foodSources.Count; i++) {
                float distToFood = Vector2.Distance(ant.position, foodSources[i].position) - foodSources[i].size;
                if (distToFood < closestFoodDistance) {
                    closestSource = foodSources[i];
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
                    return;
                }
            }
        }
    }

    private static AntData UpdateAntRotation(AntData antData) {
        AntJobData ant = antData.baseData;
        antData.sensorStrengths = new float[3] {0, 0, 0}; // Initialize sensor strengths to zero

        foreach (PheromoneData pheromone in antData.pheromonesInRange) {
            for (int i = 0; i < antData.sensorPositions.Length; i++) {
                antData.sensorStrengths[i] += CalculatePheromoneConcentration(ref antData, pheromone, i);
            }
        }
        // Normalize sensor strengths
        for (int j = 0; j < 3; j++) {
            if (antData.sensorStrengths[j] != 0) antData.sensorStrengths[j] /= antData.pheromonesInRange.Count;
        }
        bool leftSensorStrongest = (antData.sensorStrengths[1] >= antData.sensorStrengths[0] && antData.sensorStrengths[1] > antData.sensorStrengths[2]);
        bool rightSensorStrongest = (antData.sensorStrengths[2] >= antData.sensorStrengths[0] && antData.sensorStrengths[2] > antData.sensorStrengths[1]);
        // Calculate the direction of rotation based on sensor strengths
        float rotationDirection = 0f;
        if (leftSensorStrongest) {
            rotationDirection = -1;
        }
        else if (rightSensorStrongest) {
            rotationDirection = 1;
        }
        // Apply the rotation
        ant.rotation += rotationDirection * ant.rotationSpeed;
        ant.rotation = NormaliseAngle(ant.rotation);
        antData.baseData = ant;
        return antData;
    }

    private static float CalculatePheromoneConcentration (ref AntData antData, PheromoneData pheromone, int SensorIndex) {
        AntJobData ant = antData.baseData;
        float sensorDistanceToPheromone = Vector2.Distance(pheromone.position, antData.sensorPositions[SensorIndex]);
        float totalStrength = 0;
        if (sensorDistanceToPheromone <= ant.sensorRadius) {
            bool relevantPheromone = (ant.isCarryingFood == 1) ? true : (pheromone.owner.isCarryingFood == 1);
            float distanceHome = Vector2.Distance(pheromone.position, nestPosition);
            float antDistanceHome = Vector2.Distance(ant.position, nestPosition);
            if (relevantPheromone) {
                float distanceRatio = (ant.isCarryingFood == 1) ? antDistanceHome / (distanceHome + float.Epsilon) : distanceHome / (antDistanceHome + float.Epsilon);
                float directionStrength = (ant.isCarryingFood == 1) ? ((distanceHome < antDistanceHome) ? distanceRatio : -distanceRatio) :
                                                                    ((distanceHome > antDistanceHome) ? distanceRatio : -distanceRatio);
                float distanceStrength = distanceHome / directionStrength;
                float sensorDistanceStrength = sensorDistanceToPheromone / ant.sensorRadius;
                totalStrength = (distanceStrength * sensorDistanceStrength);
            }
        }
        return totalStrength;
    }

    public static float NormaliseAngle(float Angle) {
        Angle = Angle % 360f; // Normalize rotation
        if (Angle < 0) Angle += 360f; // Between 0 and 359
        return Angle;
    }

    public static void DepositPheromone(AntJobData ant) {
        int pointsToFood = ant.isCarryingFood; // Points to food if ant is carrying food
        for(int i = 0; i < pheromones.Count; i++) {
            PheromoneData pheromone = pheromones[i];
            if (Vector2.Distance(pheromone.position, ant.position) <= ant.pheromoneMergeRadius) {
                if(pointsToFood == pheromone.owner.isCarryingFood) {
                    pheromones.RemoveAt(i);
                    i--; // Adjust index after removal
                }
            }
        }
        pheromones.Add(new PheromoneData {
            position = ant.position,
            owner = ant,
            age = 0f
        });
        
    }
    
    public static FoodSource CreateNewFoodSource() {
        Vector2 pos;
        do {
            pos = new Vector2(
                UnityEngine.Random.Range(-worldSize/2, worldSize/2),
                UnityEngine.Random.Range(-worldSize/2, worldSize/2)
            );
        } while (Vector2.Distance(pos, nestPosition) < worldSize * 0.05f);
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
