#include "task.h"
#include "instance.h"
#include "../libs/vector.h"
#include "sort.h"
#include <limits>
#include <iostream>
#include <emscripten/atomic.h>

//

/* Task::Task(uint32_t id, std::function<void()> fn) :
  id(id),
  fn(fn),
  live(true),
  worldPosition{
    0,
    0,
    0
  },
  lod(0),
  priority(0)
{} */
Task::Task(uint32_t id, int priority, std::function<void()> fn) :
  id(id),
  fn(fn),
  live(true),
  worldPosition{
    0,
    0,
    0
  },
  lod(0),
  priority(priority)
{}
Task::Task(uint32_t id, const vm::vec3 &worldPosition, int lod, std::function<void()> fn) :
  id(id),
  fn(fn),
  live(true),
  worldPosition(worldPosition),
  lod(lod),
  priority(0)
{}
Task::Task(uint32_t id, const vm::vec3 &worldPosition, int lod, int priority, std::function<void()> fn) :
  id(id),
  fn(fn),
  live(true),
  worldPosition(worldPosition),
  lod(lod),
  priority(priority)
{}

Task::~Task() {}

void Task::run() {
  fn();
}
void Task::cancel() {
  live.store(false);
}

int Task::getPriority() const {
  return priority;
}
/* Sphere Task::getSphere() const {
  return Sphere{
    Vec{
      worldPosition.x,
      worldPosition.y,
      worldPosition.z
    },
    std::sqrt(
      halfSize.x * halfSize.x +
      halfSize.y * halfSize.y +
      halfSize.z * halfSize.z
    )
  };
} */
Box3 Task::getBox() const {
  return Box3{
    Vec{
      worldPosition.x,
      worldPosition.y,
      worldPosition.z
    },
    Vec{
      worldPosition.x + lod,
      worldPosition.y + lod,
      worldPosition.z + lod
    }
  };
}

//

TaskQueue::TaskQueue() {}
TaskQueue::~TaskQueue() {
  // EM_ASM({
  //   console.log('task queue destructor');
  // });
  abort();
}

void TaskQueue::pushTask(Task *task) {
  Frustum frustum = getFrustum();
  const float taskDistance = getTaskDistance(task, frustum);
  {
    std::unique_lock<Mutex> lock(taskMutex);

    bool found = false;
    for (size_t i = 0; i < tasks.size(); i++) {
      auto iter = tasks.begin() + i;
      const float taskDistance2 = getTaskDistance(*iter, frustum);
      if (taskDistance < taskDistance2) {
        tasks.insert(iter, task);
        found = true;
        break;
      }
    }
    if (!found) {
      tasks.push_back(task);
    }
  }
  taskSemaphore.signal();
}

Task *TaskQueue::popLockTask() {
  taskSemaphore.wait();
  
  Task *task;
  {
    std::unique_lock<Mutex> lock(taskMutex);

    task = tasks.front();
    tasks.pop_front();
  }
  
  return task;
}
void TaskQueue::cancelTask(uint32_t taskId) {
  std::unique_lock<Mutex> lock(taskMutex);
  for (auto it = tasks.begin(); it != tasks.end(); it++) {
    Task *task = (*it);
    if (task->id == taskId && task->live) {
      task->cancel();
      break;
    }
  }
}
void TaskQueue::runLoop() {
    for (;;) {
      Task *task = popLockTask();
      if (!task) {
        std::cout << "failed to pop task" << std::endl;
        abort();
      }
      if (task->live) {
        task->run();
      }
      delete task;
    }

  std::cout << "main loop exited" << std::endl;
  abort();
}

void TaskQueue::setCamera(const vm::vec3 &worldPosition, const vm::vec3 &cameraPosition, const Quat &cameraQuaternion, const std::array<float, 16> &projectionMatrix) {
  std::unique_lock<Mutex> lock(taskMutex);

  this->worldPosition = worldPosition;
  this->cameraPosition = cameraPosition;
  this->cameraQuaternion = cameraQuaternion;
  this->projectionMatrix = projectionMatrix;

  sortTasksInternal();
}
Frustum TaskQueue::getFrustum() {
  Matrix matrixWorld(
    Vec{
      cameraPosition.x,
      cameraPosition.y,
      cameraPosition.z
    },
    Quat{
      cameraQuaternion.x,
      cameraQuaternion.y,
      cameraQuaternion.z,
      cameraQuaternion.w
    },
    Vec{1, 1, 1}
  );
  Matrix matrixWorldInverse(matrixWorld);
  matrixWorldInverse.invert();
  Frustum frustum = Frustum::fromMatrix(
    Matrix::fromArray(projectionMatrix.data()) *= matrixWorldInverse
  );
  return frustum;
}
float TaskQueue::getTaskDistance(Task *task, const Frustum &frustum) {
  double distance = vm::length(task->worldPosition - worldPosition);

  const float halfLod = (float)task->lod / 2.0f;
  Box3 box(
    Vec{
      task->worldPosition.x - halfLod,
      task->worldPosition.y - halfLod,
      task->worldPosition.z - halfLod
    },
    Vec{
      task->worldPosition.x + halfLod,
      task->worldPosition.y + halfLod,
      task->worldPosition.z + halfLod
    }
  );
  if (!frustum.intersectsBox(box)) {
    distance += frustumCullDistancePenalty;
  }

  distance += task->priority * priorityDistancePenalty;
  return distance;
}
void TaskQueue::sortTasksInternal() {
  const vm::vec3 &worldPosition = this->worldPosition;
  Frustum frustum = getFrustum();

  sort<Task *>(tasks, worldPosition, frustum);
}