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
        AntSimulationCore.Initialize(colonySize, nestPosition, 3, worldSize, 0.5);
        if(colonySize > 256) CreateAntsWithJobSystem(colonySize);
        else AntSimulationCore.CreateAnts(colonySize);
        SwitchRenderer(currentRenderMode);
    }

    void Update() {
        if (Input.GetMouseButtonDown(0)) {
            Vector2 mousePos = Camera.main.ScreenToWorldPoint(Input.mousePosition);
            AntSimulationCore.foodSources.Add(new FoodSource {
                position = mousePos,
                mass = UnityEngine.Random.Range(50f, 200f),
                size = UnityEngine.Random.Range(0.5f, 2f)
            });
        }
        if(AntSimulationCore.pheromones.Count > 256) UpdatePheromonesWithJobSystem(Time.deltaTime);
        else AntSimulationCore.UpdatePheromoneAges(Time.deltaTime);
        AntSimulationCore.Update(Time.deltaTime);
        UpdateSimulationData();
        
        if(currentRenderMode != RenderMode.Gizmos) {
            _activeRenderer.Draw(_simulationData);
        }
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
        using (var removeQueue = new NativeQueue<int>(Allocator.TempJob)) {
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
            // Create a random number generator with a unique seed per ant
            var random = new Unity.Mathematics.Random(Seed + (uint)index);

            // Generate random position within nest radius
            float2 randomOffset = random.NextFloat2(-NestSize / 2, NestSize / 2);
            float2 position = NestPosition + randomOffset;

            // Generate random rotation
            float rotation = random.NextFloat(0, 360f);

            // Generate random values for ant properties using Unity.Mathematics.Random
            float PheromoneIntervalValue = random.NextFloat(0f, 1f);
            float PheromoneDecayRate = random.NextFloat(0f, 0.1f);
            float PheromoneMaxAge = (1 - PheromoneIntervalValue) / PheromoneDecayRate;
            float PheromoneMergeRadius = random.NextFloat(0f, 1f) / 10;
            float Speed = random.NextFloat(0f, 1f) * 2;
            float RotationSpeed = random.NextFloat(0f, 360f) / 2;
            float MaxEnergy = (PheromoneDecayRate + PheromoneMergeRadius + Speed + (RotationSpeed / 360)) / (PheromoneIntervalValue + PheromoneMaxAge) * 150;

            // Create the ant data
            AntJobData ant = new AntJobData {
                position = position,
                rotation = rotation,
                isCarryingFood = 0, // false
                pheromoneTimer = 0f,
                pheromoneIntervalValue = PheromoneIntervalValue,
                pheromoneDecayRate = PheromoneDecayRate,
                pheromoneMaxAge = PheromoneMaxAge,
                pheromoneMergeRadius = PheromoneMergeRadius,
                sensorRadius = random.NextFloat(0f, 1f) * 2,
                smellRadius = random.NextFloat(1f, 2f) * 2,
                speed = Speed,
                rotationSpeed = RotationSpeed,
                maxEnergy = MaxEnergy,
                energy = MaxEnergy,
                energyConsumptionRate = (PheromoneMaxAge * Speed * (RotationSpeed / 180)) / 10
            };

            // Store the ant
            Ants[index] = ant;
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
            if (pheromone.age >= pheromone.owner.pheromoneMaxAge) RemoveQueueWriter.Enqueue(index); // Now using the properly defined field
            PheromoneData[index] = pheromone;
        }
    }
}

[System.Serializable]
public struct AntSimulationData {
    public float SimulationTime;
    public float2 NestPosition;
    public double NestSize;
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
    public float pheromoneMergeRadius; // Radius to check for nearby pheromones

    public float sensorRadius;
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
    public double mass;
    public double size;
}

public static class AntSimulationCore {
    public static float simulationTime = 0f;
    public static float2 nestPosition;
    public static double nestSize;
    public static float worldSize = 50f;
    public static List<AntData> ants = new List<AntData>();
    public static List<PheromoneData> pheromones = new List<PheromoneData>();
    public static List<FoodSource> foodSources = new List<FoodSource>();

    public static void Initialize(int colonySize, Vector2 nestPos, int foodSourceCount, float WorldSize, double NestSize) {
        nestPosition = nestPos;
        ants = new List<AntData>();
        pheromones = new List<PheromoneData>();
        foodSources = new List<FoodSource>();
        simulationTime = 0f;
        worldSize = WorldSize;
        nestSize = NestSize;
        
        // Create food sources
        for(int i = 0; i < foodSourceCount; i++) {
            Vector2 pos;
            do {
                pos = new Vector2(
                    UnityEngine.Random.Range(-worldSize/2, worldSize/2),
                    UnityEngine.Random.Range(-worldSize/2, worldSize/2)
                );
            } while (Vector2.Distance(pos, nestPos) < worldSize * 0.01); // Ensure food isn't right at nest
            float massa = UnityEngine.Random.Range(500f, 20000f);
            foodSources.Add(new FoodSource {
                position = pos,
                mass = massa,
                size = Mathf.Max(massa / 1000 , 0.01f)
            });
        }
    }

    public static void CreateAnts(int ColonySize) {
        for(int i = 0; i < ColonySize; i++) {
            float PheromoneIntervalValue = UnityEngine.Random.Range(0f, 1f);
            float PheromoneDecayRate = UnityEngine.Random.Range(0f, 0.1f);
            float PheromoneMaxAge = (1 - PheromoneIntervalValue) / PheromoneDecayRate;
            float PheromoneMergeRadius = UnityEngine.Random.Range(0f, 1f) / 10;
            float Speed = UnityEngine.Random.Range(0f, 1f) * 2;
            float RotationSpeed = UnityEngine.Random.Range(0f, 360f) / 2;
            float MaxEnergy = (PheromoneDecayRate + PheromoneMergeRadius + Speed + (RotationSpeed / 360)) / (PheromoneIntervalValue + PheromoneMaxAge) * 150;
            ants.Add(new AntData {
                baseData = new AntJobData {
                    position = new Vector2(nestPosition.x + UnityEngine.Random.Range((float)-nestSize / 2, (float)nestSize / 2), nestPosition.y + UnityEngine.Random.Range((float)-nestSize / 2, (float)nestSize / 2)),
                    rotation = UnityEngine.Random.Range(0f, 360f),
                    isCarryingFood = 0,

                    pheromoneTimer = 0,
                    pheromoneIntervalValue = PheromoneIntervalValue,
                    pheromoneDecayRate = PheromoneDecayRate,
                    pheromoneMaxAge = PheromoneMaxAge,
                    pheromoneMergeRadius = PheromoneMergeRadius,

                    sensorRadius = UnityEngine.Random.Range(0f, 1f) * 2,
                    smellRadius = UnityEngine.Random.Range(1f, 2f) * 2,

                    speed = Speed,
                    rotationSpeed = RotationSpeed,
                    maxEnergy = MaxEnergy,
                    energy = MaxEnergy,
                    energyConsumptionRate = (PheromoneMaxAge * Speed * (RotationSpeed / 180)) / 10
                },
                pheromonesInRange = new List<PheromoneData>()
            });
        }
    }

    public static void UpdatePheromoneAges(float deltaTime) {
        if(pheromones.Count > 0) {
            for (int i = pheromones.Count - 1; i >= 0; i--) {
                PheromoneData pheromone = pheromones[i];
                pheromone.age += pheromone.owner.pheromoneDecayRate * deltaTime;
                if (pheromone.age >= pheromone.owner.pheromoneMaxAge) {
                    pheromones.RemoveAt(i); // Remove expired pheromones
                } else {
                    pheromones[i] = pheromone;
                }
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
            else if (Vector2.Distance(ant.position, nestPosition) <= nestSize * 2) {
                // Ant delivered food!
                ant.isCarryingFood = 0;
                ant.energy = ant.maxEnergy; // Reward energy
                
                // Turn around to leave nest
                ant.rotation += 180f;
            }
            // Update ant direction based on pheromones
            UpdateAntDirection(ref ant, deltaTime);
            // Move ant
            ant.position += new float2(Mathf.Cos(ant.rotation * Mathf.Deg2Rad), Mathf.Sin(ant.rotation * Mathf.Deg2Rad)) * ant.speed * deltaTime;
            if (ant.position.x < -worldSize/2 || ant.position.x > worldSize/2 || ant.position.y < -worldSize/2 || ant.position.y > worldSize/2) {
                ant.position = nestPosition;
                ant.rotation += 180f;
            }
            antData.baseData = ant;
            ants[i] = antData;
        }
    }

    public static void Update(float deltaTime) {
        simulationTime += deltaTime;
        UpdateAnts(deltaTime);
        // Check for depleted food sources
        if (foodSources.Count > 0) {
            for (int i = foodSources.Count - 1; i >= 0; i--) {
                if (foodSources[i].mass <= 0) {
                    foodSources.RemoveAt(i);
                }
            }
        }
    }
    
    private static void UpdateAntDirection(ref AntJobData ant, float deltaTime) {
        // Special case: If carrying food and near nest, steer toward nest
        if (ant.isCarryingFood == 1) {
            float distToNest = Vector2.Distance(ant.position, nestPosition);
            if (distToNest < nestSize * 1.75) { 
                float2 dirToNest = (float2)nestPosition - ant.position;
                float targetAngle = Mathf.Atan2(dirToNest.y, dirToNest.x) * Mathf.Rad2Deg;
                float angleDiff = Mathf.DeltaAngle(ant.rotation, targetAngle);
                ant.rotation += Mathf.Clamp(angleDiff, -ant.rotationSpeed, ant.rotationSpeed);
                return;
            }
        }
        // Special case: If not carrying food and near food, steer toward food
        if (ant.isCarryingFood == 0) {
            float closestFoodDistance = float.MaxValue;
            FoodSource? closestSource = null;
            for(int i = 0; i < foodSources.Count; i++){
                float distToFood = Vector2.Distance(ant.position, foodSources[i].position);
                if(distToFood < closestFoodDistance) {
                    closestSource = foodSources[i];
                    closestFoodDistance = distToFood;
                }
            }
            if(closestSource != null) {
                if (closestFoodDistance <= ant.smellRadius + closestSource.Value.size) {
                    Vector2 dirToFood = (float2)closestSource.Value.position - ant.position;
                    float targetAngle = Mathf.Atan2(dirToFood.y, dirToFood.x) * Mathf.Rad2Deg;
                    float angleDiff = Mathf.DeltaAngle(ant.rotation, targetAngle);
                    ant.rotation += Mathf.Clamp(angleDiff, -ant.rotationSpeed, ant.rotationSpeed);
                    return;
                }
            }
        }
        // Create sensors at different angles around the ant
        Vector2[] sensorPositions = new Vector2[3];
        float[] sensorStrengths = new float[3];
        // Front, left and right sensors
        float sensorAngle = 60f;
        float antDirection = ant.rotation * Mathf.Deg2Rad;
        sensorPositions[0] = ant.position + new float2(
            Mathf.Cos(antDirection) * ant.sensorRadius,
            Mathf.Sin(antDirection) * ant.sensorRadius
        );
        sensorPositions[1] = ant.position + new float2(
            Mathf.Cos(antDirection - sensorAngle * Mathf.Deg2Rad) * ant.sensorRadius,
            Mathf.Sin(antDirection - sensorAngle * Mathf.Deg2Rad) * ant.sensorRadius
        );
        sensorPositions[2] = ant.position + new float2(
            Mathf.Cos(antDirection + sensorAngle * Mathf.Deg2Rad) * ant.sensorRadius,
            Mathf.Sin(antDirection + sensorAngle * Mathf.Deg2Rad) * ant.sensorRadius
        );
        
        // Calculate pheromone concentrations at each sensor
        for (int i = 0; i < sensorPositions.Length; i++) {
            sensorStrengths[i] = 0f;
            foreach (var pheromone in pheromones) {
                // Check if this is the right type of pheromone for the ant's state
                bool relevantPheromone = (ant.isCarryingFood == 1) ? true : (pheromone.owner.isCarryingFood == 1);
                float distance = Vector2.Distance(sensorPositions[i], pheromone.position);
                float distanceHome = Vector2.Distance(pheromone.position, nestPosition);
                float antDistanceHome = Vector2.Distance(ant.position, nestPosition);
                float directionStrength = 0;
                if (relevantPheromone && ant.isCarryingFood == 0) {
                    if (distance <= ant.sensorRadius) {
                        directionStrength = (distanceHome > antDistanceHome) ? distanceHome / antDistanceHome : -distanceHome / antDistanceHome;
                        sensorStrengths[i] += directionStrength;
                    }
                }
                if (relevantPheromone && ant.isCarryingFood == 1) {
                    if (distance <= ant.sensorRadius) {
                        directionStrength = (distanceHome < antDistanceHome) ? antDistanceHome / distanceHome : -antDistanceHome / distanceHome;
                        sensorStrengths[i] += directionStrength;
                    }
                }
            }
        }
        
        // Check if any sensor detected pheromones
        bool pheromoneDetected = sensorStrengths[0] > 0 || sensorStrengths[1] > 0 || sensorStrengths[2] > 0;
        if (pheromoneDetected) {
            // Determine turning direction based on sensor readings
            if (sensorStrengths[1] > sensorStrengths[0] && sensorStrengths[1] > sensorStrengths[2]) {
                // Left sensor strongest, turn left
                ant.rotation -= ant.rotationSpeed;
            } 
            else if (sensorStrengths[2] > sensorStrengths[0] && sensorStrengths[2] > sensorStrengths[1]) {
                // Right sensor strongest, turn right
                ant.rotation += ant.rotationSpeed;
            }
            // Otherwise continue straight (center sensor strongest)
            return;
        }
        
        // No pheromones detected, use random movement with occasional direction changes
        if (UnityEngine.Random.value < 0.02f) { // 2% chance per frame to change direction
            ant.rotation += UnityEngine.Random.Range(-ant.rotationSpeed, ant.rotationSpeed);
        }
    }

    public static void DepositPheromone(AntJobData ant) {
        bool pointsToFood = ant.isCarryingFood ==1; // Points to food if ant is carrying food
        bool shouldDrop = true;

        // First loop: Remove competing pheromones
        for(int i = 0; i < pheromones.Count; i++) {
            PheromoneData pheromone = pheromones[i];
            if (Vector2.Distance(pheromone.position, ant.position) < ant.pheromoneMergeRadius) {
                if(pointsToFood || pheromone.owner.isCarryingFood == 0) {
                    pheromones.RemoveAt(i);
                    i--; // Adjust index after removal
                }
                else if((pheromone.owner.isCarryingFood == 1) == pointsToFood) {
                    pheromone.age = 0;
                    pheromone.position = ant.position;
                    pheromones[i] = pheromone;
                }
            }
        }
        
        // Second loop: Check if we should drop a new pheromone
        for(int j = 0; j < pheromones.Count; j++) {
            PheromoneData pheromone = pheromones[j];
            if (Vector2.Distance(pheromone.position, ant.position) < ant.pheromoneMergeRadius / 2 && !pointsToFood) {
                shouldDrop = false;
                break;
            }
        }
        
        // Add new pheromone if appropriate
        if(shouldDrop) {
            pheromones.Add(new PheromoneData {
                position = ant.position,
                owner = ant,
                age = 0f
            });
        }
        
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
        foreach(var ant in ants) {
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
