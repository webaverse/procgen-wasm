#include <emscripten.h>
// #include "DualContouring/tracker.h"
// #include "DualContouring/main.h"
#include "procgen.h"
#include "meshoptimizer/meshoptimizer.h"

extern "C" {

EMSCRIPTEN_KEEPALIVE void initialize() {
    ProcGen::initialize();
}

// 

EMSCRIPTEN_KEEPALIVE PGInstance *createInstance(int seed, int chunkSize) {
    // std::cout << "create instance " << seed << " " << chunkSize << std::endl;
    return ProcGen::createInstance(seed, chunkSize);
}
EMSCRIPTEN_KEEPALIVE void destroyInstance(PGInstance *instance) {
    ProcGen::destroyInstance(instance);
}

//

/* EMSCRIPTEN_KEEPALIVE void getHeightfieldRangeAsync(DCInstance *inst, uint32_t taskId, int x, int z, int w, int h, int lod, float *heights, int priority) {
    vm::ivec2 worldPositionXZ{x, z};
    vm::ivec2 sizeXZ{w, h};
    inst->getHeightfieldRangeAsync(taskId, worldPositionXZ, sizeXZ, lod, heights, priority);
}
EMSCRIPTEN_KEEPALIVE void getLightRangeAsync(DCInstance *inst, uint32_t taskId, int x, int y, int z, int w, int h, int d, int lod, uint8_t *skylights, uint8_t *aos, int priority) {
    vm::ivec3 worldPosition{x, y, z};
    vm::ivec3 size{w, h, d};
    inst->getLightRangeAsync(taskId, worldPosition, size, lod, skylights, aos, priority);
} */

// 

/* EMSCRIPTEN_KEEPALIVE void createGrassSplatAsync(DCInstance *inst, uint32_t taskId, int x, int z, int lod, int priority) {
    inst->createGrassSplatAsync(taskId, vm::ivec2{x, z}, lod, priority);
}
EMSCRIPTEN_KEEPALIVE void createVegetationSplatAsync(DCInstance *inst, uint32_t taskId, int x, int z, int lod, int priority) {
    inst->createVegetationSplatAsync(taskId, vm::ivec2{x, z}, lod, priority);
}
EMSCRIPTEN_KEEPALIVE void createMobSplatAsync(DCInstance *inst, uint32_t taskId, int x, int z, int lod, int priority) {
    inst->createMobSplatAsync(taskId, vm::ivec2{x, z}, lod, priority);
} */

//

EMSCRIPTEN_KEEPALIVE void createChunkMeshAsync(
    PGInstance *inst,
    uint32_t taskId,
    int x, int z,
    int lod,
    int *lodArray,
    int generateFlags,
    int numVegetationInstances,
    int numRockInstances,
    int numGrassInstances,
    int numPoiInstances
) {
    vm::ivec2 worldPosition{x, z};
    std::array<int, 2> lodArray2{
        lodArray[0],
        lodArray[1]
    };
    inst->createChunkMeshAsync(
        taskId,
        worldPosition,
        lod,
        lodArray2,
        generateFlags,
        numVegetationInstances,
        numRockInstances,
        numGrassInstances,
        numPoiInstances
    );
}
EMSCRIPTEN_KEEPALIVE void createBarrierMeshAsync(
    PGInstance *inst,
    uint32_t taskId,
    int x, int z,
    int minLod,
    int maxLod
) {
    vm::ivec2 worldPosition{x, z};
    inst->createBarrierMeshAsync(
        taskId,
        worldPosition,
        minLod,
        maxLod
    );
}
/* EMSCRIPTEN_KEEPALIVE void createChunkGrassAsync(PGInstance *inst, uint32_t taskId, int x, int z, int lod, int numGrassInstances) {
    inst->createChunkGrassAsync(taskId, vm::ivec2{x, z}, lod, numGrassInstances);
} */
/* EMSCRIPTEN_KEEPALIVE void createChunkVegetationAsync(PGInstance *inst, uint32_t taskId, int x, int z, int lod, int numVegetationInstances) {
    inst->createChunkVegetationAsync(taskId, vm::ivec2{x, z}, lod, numVegetationInstances);
} */

/* EMSCRIPTEN_KEEPALIVE void createLiquidChunkMeshAsync(DCInstance *inst, uint32_t taskId, int x, int y, int z, int *lodArray) {
    inst->createLiquidChunkMeshAsync(taskId, vm::ivec3{x, y, z}, lodArray);
} */

//

/* EMSCRIPTEN_KEEPALIVE bool drawSphereDamage(DCInstance *inst, float x, float y, float z, float radius, float *outPositions, unsigned int *outPositionsCount) {
    return inst->drawSphereDamage(x, y, z, radius, outPositions, outPositionsCount, 1);
}

EMSCRIPTEN_KEEPALIVE bool eraseSphereDamage(DCInstance *inst, float x, float y, float z, float radius, float *outPositions, unsigned int *outPositionsCount, float *outDamages) {
    return inst->eraseSphereDamage(x, y, z, radius, outPositions, outPositionsCount, outDamages, 1);
}

EMSCRIPTEN_KEEPALIVE bool drawCubeDamage(
    DCInstance *inst,
    float x, float y, float z,
    float qx, float qy, float qz, float qw,
    float sx, float sy, float sz,
    float *outPositions,
    unsigned int *outPositionsCount,
    float *outDamages
) {
    return inst->drawCubeDamage(
        x, y, z,
        qx, qy, qz, qw,
        sx, sy, sz,
        outPositions,
        outPositionsCount,
        outDamages,
        1
    );
}

EMSCRIPTEN_KEEPALIVE bool eraseCubeDamage(
    DCInstance *inst,
    float x, float y, float z,
    float qx, float qy, float qz, float qw,
    float sx, float sy, float sz,
    float *outPositions,
    unsigned int *outPositionsCount,
    float *outDamages
) {
    return inst->eraseCubeDamage(
        x, y, z,
        qx, qy, qz, qw,
        sx, sy, sz,
        outPositions,
        outPositionsCount,
        outDamages,
        1
    );
} */

/* EMSCRIPTEN_KEEPALIVE void injectDamage(DCInstance *inst, float x, float y, float z, float *damageBuffer) {
    inst->injectDamage(x, y, z, damageBuffer, 1);
} */

//

EMSCRIPTEN_KEEPALIVE void cancelTask(PGInstance *inst, uint32_t taskId) {
    ProcGen::taskQueue.cancelTask(taskId);
    ProcGen::resultQueue.cancelPromise(taskId);
}

//

EMSCRIPTEN_KEEPALIVE void setCamera(PGInstance *inst, float *worldPosition, float *cameraPosition, float *cameraQuaternion, float *projectionMatrix) {
    vm::vec3 _worldPosition{
        worldPosition[0],
        worldPosition[1],
        worldPosition[2]
    };
    vm::vec3 _cameraPosition{
        cameraPosition[0],
        cameraPosition[1],
        cameraPosition[2]
    };
    Quat _cameraQuaternion{
        cameraQuaternion[0],
        cameraQuaternion[1],
        cameraQuaternion[2],
        cameraQuaternion[3]
    };
    std::array<float, 16> _projectionMatrix;
    memcpy(&_projectionMatrix[0], projectionMatrix, sizeof(_projectionMatrix));
    inst->setCamera(
        _worldPosition,
        _cameraPosition,
        _cameraQuaternion,
        _projectionMatrix
    );
}

//

EMSCRIPTEN_KEEPALIVE Tracker *createTracker(
    PGInstance *inst
) {
    // std::cout << "create tracker 1 " << (void *)inst << " " << lods << " " << lod1Range << std::endl;
    Tracker *tracker = new Tracker(inst);
    // std::cout << "create tracker 2 " << (void *)inst << " " << lods << " " << lod1Range << std::endl;
    return tracker;
}

EMSCRIPTEN_KEEPALIVE void trackerUpdateAsync(
    PGInstance *inst,
    uint32_t taskId,
    Tracker *tracker,
    float *position,
    int minLod,
    int maxLod,
    int lod1Range,
    int priority
) {
    vm::vec3 worldPosition{
        position[0],
        position[1],
        position[2],
    };
    inst->trackerUpdateAsync(taskId, tracker, worldPosition, minLod, maxLod, lod1Range, priority);
}

EMSCRIPTEN_KEEPALIVE void destroyTracker(PGInstance *inst, Tracker *tracker) {
    delete tracker;
}

//

EMSCRIPTEN_KEEPALIVE void *doMalloc(size_t size, PGInstance *inst) {
    void* ptr = malloc(size);
    return ptr;
}

EMSCRIPTEN_KEEPALIVE void doFree(void *ptr, PGInstance *inst) {
    free(ptr);
}

//

EMSCRIPTEN_KEEPALIVE void runLoop() {
    ProcGen::runLoop();
}

//

EMSCRIPTEN_KEEPALIVE void meshoptSimplify(
  const unsigned int* indices,
  size_t index_count,
  const float* vertex_positions,
  size_t vertex_count,
  // const float* uvs,
  size_t target_index_count,
  float target_error,
  unsigned int **destination,
  unsigned int *numDestinations
) {
  *destination = (unsigned int *)malloc(sizeof(unsigned int) * index_count);
  
  const size_t vertex_positions_stride = 3 * sizeof(float);
  constexpr unsigned int options = 0;
  float result_error;
  *numDestinations = meshopt_simplify(*destination, indices, index_count, vertex_positions, vertex_count, vertex_positions_stride, target_index_count, target_error, options, &result_error);
}

//

int main() {
    ProcGen::start();
    return 0;
}

} // extern "C"