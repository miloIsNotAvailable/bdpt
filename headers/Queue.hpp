#ifndef QUEUE_H
#define QUEUE_H

#include <iostream>
#include <thread>
#include <functional>
#include <condition_variable>
#include "tiny_obj_loader.h"

class Queue
{
private:
    /* data */
public:
    Queue(/* args */);
    ~Queue();

    void run();
    void add( int *a );
};

class Counter {
    
    private:
        int count;
    
    public:
        int id;
    
        Counter( int id ) {
            this->id=id;
            this->count=0;
        }

        void increment() {
            this->count++;
        }
        int get() {
            return this->count;
        }

};

template <typename T> class Thread {
    public: 
        std::thread::id id;
        T *data;
        std::thread thread;
        std::condition_variable cv;

        Thread(int threadId, T* inputData, std::function<void(T*, int)> func)
        : id(threadId), data(inputData) {
            this->thread = std::thread(func, data);
            this->id = this->thread.get_id();
            this->cv;
        }

        Thread(const Thread&) = delete;
        Thread& operator=(const Thread&) = delete;
        Thread(Thread&&) = default;
        Thread& operator=(Thread&&) = default;

        void join() {
            if (!this->thread.joinable()) throw std::runtime_error("error");
            this->thread.join();
        }
};

typedef struct {
    std::vector<float> vertices;
    std::vector<float> normals;
    std::vector<float> texcoords;
    std::vector<int> v_indices;
    std::vector<int> vn_indices;
    std::vector<int> vt_indices;
  
    std::vector<tinyobj::material_t> materials;
  
  } MyMesh;  


#endif