#ifndef PROMISE_H
#define PROMISE_H

#include "../libs/vectorMath.h"
#include "sync.h"
#include "result.h"

#include <vector>
#include <deque>
#include <atomic>
#include <emscripten.h>

//

class Promise {
public:
    uint32_t id;
    ResultQueue *resultQueue;

    std::atomic<bool> live;

    Promise(uint32_t id, ResultQueue *resultQueue);
    ~Promise();

    bool resolve(void *value = nullptr);
    bool kill();
};

#endif // PROMISE_H
