#include "procgen.h"
#include "instance.h"
#include <pthread.h>
#include <thread>

namespace ProcGen
{
    TaskQueue taskQueue;
    ResultQueue resultQueue;

    pthread_t parentThreadId;

    void start() {
        // threads.reserve(numThreads);
        /* EM_ASM({
            console.log('main thread', $0, $1);
        }, pthread_self(), numThreads); */
        for (int i = 0; i < numThreads; i++) {
            // std::cout << "create thread" << std::endl;
            std::thread([]() -> void {
                /* EM_ASM({
                    console.log('worker thread', $0);
                }, pthread_self()); */
                runLoop();
            }).detach();
        }
    }

    void runLoop() {
        taskQueue.runLoop();
    }

    void initialize()
    {
        parentThreadId = pthread_self();
    }

    PGInstance *createInstance(int seed, int chunkSize) {
        PGInstance *instance = new PGInstance(seed, chunkSize);
        return instance;
    }
    void destroyInstance(PGInstance *instance) {
        delete instance;
    }
}