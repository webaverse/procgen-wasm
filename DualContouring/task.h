#ifndef TASK_H
#define TASK_H

#include "vectorMath.h"
#include "sync.h"
#include "lock.h"
#include <vector>
#include <deque>
// #include <semaphore>
#include <atomic>

//

class DCInstance;

//

class Task {
public:
    MultiChunkLock multiChunkLock;
    std::function<void()> fn;

    Task(MultiChunkLock &&multiChunkLock, std::function<void()> &&fn);
    ~Task();

    bool tryLock();
    void lock();
    void unlock();
    void run();
    std::pair<bool, void *> tryLockRun();
};

//

class TaskQueue {
public:
    DCInstance *inst;
    std::deque<Task *> tasks;
    // std::atomic<size_t> numTasks;
    Mutex taskMutex;
    Semaphore taskSemaphore;

    TaskQueue();
    ~TaskQueue();
    
    void pushTask(Task *task);
    Task *popLockTask();
    void runLoop();
};

#endif // TASK_H
