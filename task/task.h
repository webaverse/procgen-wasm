#ifndef TASK_H
#define TASK_H

#include "../libs/vectorMath.h"
#include "../libs/vector.h"
#include "sync.h"

#include <array>
#include <vector>
#include <deque>
#include <atomic>
#include <emscripten.h>


class Task {
public:
    uint32_t id;
    std::function<void()> fn;
    std::atomic<bool> live;

    vm::vec3 worldPosition;
    int lod;
    int priority;

    Task(uint32_t id, int priority, std::function<void()> fn);
    Task(uint32_t id, const vm::vec3 &worldPosition, int lod, std::function<void()> fn);
    Task(uint32_t id, const vm::vec3 &worldPosition, int lod, int priority, std::function<void()> fn);
    ~Task();

    void run();
    void cancel();

    int getPriority() const;
    Box3 getBox() const;
};

//

class TaskQueue {
public:
    // DCInstance *inst;
    std::deque<Task *> tasks;
    
    Mutex taskMutex;
    Semaphore taskSemaphore;

    vm::vec3 worldPosition;
    vm::vec3 cameraPosition;
    Quat cameraQuaternion;
    std::array<float, 16> projectionMatrix;

    TaskQueue();
    ~TaskQueue();
    
    void pushTask(Task *task);
    Task *popLockTask();
    void runLoop();
    void cancelTask(uint32_t taskId);

    void setCamera(const vm::vec3 &worldPosition, const vm::vec3 &cameraPosition, const Quat &cameraQuaternion, const std::array<float, 16> &projectionMatrix);
    Frustum getFrustum();
    float getTaskDistance(Task *task, const Frustum &frustum);

    void sortTasksInternal();
};

#endif // TASK_H
