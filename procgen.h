#ifndef PROCGEN_H
#define PROCGEN_H

#include <iostream>
#include <vector>
#include <cstdint>
#include <ctime>
#include <string.h>
#include <memory>
#include "instance.h"
#include "noises.h"
#include "tracker.h"
#include "result.h"
#include "vectorMath.h"
#include "constants.h"

//

class Noises;
class PGInstance;

//

namespace ProcGen {
    // globals
    extern Noises *noises;
    extern TaskQueue taskQueue;
    extern ResultQueue resultQueue;

    extern pthread_t parentThreadId;

    // initialization
    void initialize(int newChunkSize, int seed);
    
    // instances
    PGInstance *createInstance(int seed, int chunkSize);
    void destroyInstance(PGInstance *instance);

    float getComputedBiomeHeight(unsigned char b, const vm::vec2 &worldPosition);

    // threads
    void start();
    void runLoop();
};

#endif // PROCGEN_H
