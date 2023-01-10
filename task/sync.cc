#include "sync.h"
#include <iostream>
#include <emscripten.h>

//

Mutex::Mutex() {}
Mutex::Mutex(bool locked) {
  if (locked) {
    mutex.try_lock();
  }
}
Mutex::~Mutex() {}
void Mutex::lock() {
  mutex.lock();
}
void Mutex::unlock() {
  mutex.unlock();
}
bool Mutex::try_lock() {
  return mutex.try_lock();
}

Semaphore::Semaphore(int count) : count(count) {}
void Semaphore::signal() {
    std::unique_lock<std::mutex> lock(mtx);
    count++;
    cv.notify_one();
}
void Semaphore::wait() {
    std::unique_lock<std::mutex> lock(mtx);
    while (count == 0) {
        cv.wait(lock);
    }
    count--;
}