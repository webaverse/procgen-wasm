#ifndef PROCGEN_H
#define PROCGEN_H

#include <iostream>
#include <vector>
#include <cstdint>
#include <ctime>
#include <string.h>
#include <memory>

#include "instance.h"
#include "generation/noises.h"
#include "task/tracker.h"
#include "task/task.h"
#include "task/result.h"
#include "libs/vectorMath.h"
#include "constants.h"

//

class Noises;
class PGInstance;

//


namespace ProcGen {
    // globals
    extern TaskQueue taskQueue;
    extern ResultQueue resultQueue;

    extern pthread_t parentThreadId;

    // initialization
    void initialize();
    
    // instances
    PGInstance *createInstance(int seed, int chunkSize);
    void destroyInstance(PGInstance *instance);

    // threads
    void start();
    void runLoop();
};

#endif // PROCGEN_H
