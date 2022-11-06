#ifndef __MEM_H__
#define __MEM_H__

#include <cstdlib>
#include <unordered_map>
#include <string>
#include <sstream>
#include <emscripten.h>

struct AllocationData {
    size_t size;
    std::string location;

    std::string getDataString() {
        std::string dataString = "";

        dataString += "( ";
        dataString += location;
        dataString += ", ";
        dataString += std::to_string(size);
        dataString += " )";

        return dataString;
    }
};

struct MemoryManager {
    std::unordered_map<void*, AllocationData> allocationMap;

    MemoryManager() {
        // EM_ASM(
        //     console.log("MemoryManager Init");
        // );
    }

    ~MemoryManager() {
        // EM_ASM(
        //     console.log("MemoryManager Destroy");
        // );
    }

    template<typename T>
    void allocate(T ptr, size_t size, const std::string &str) {
        // if(ptr) {
        //     AllocationData a{size, str};
        //     void *vptr = (void*)ptr;
        //     allocationMap.insert(std::make_pair(vptr, a));
        // }
    }

    template<typename T>
    void free(T ptr) {
        // if(ptr) {
        //     void *vptr = (void*)ptr;
        //     allocationMap.erase(vptr);
        // }
    }

    void getLeaked() {
        // size_t leaked = 0;
        // std::string str = "Leaked memory locations: ";

        // for (auto &pair : allocationMap) {
        //     AllocationData &a = pair.second;
        //     leaked += a.size;
        //     str += a.getDataString() + ", ";
        // }

        // EM_ASM({
        //     console.log('Leaked Size: ', $0);
        // }, leaked);

        // const char* ccstr = str.c_str();
        // EM_ASM({
        //     console.log(UTF8ToString($0));
        // }, ccstr);
    }
};

#endif // __MEM_H__