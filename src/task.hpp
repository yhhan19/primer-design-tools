#pragma once
#include "utility.hpp"

struct Task {
    std::string contig;
    std::string seq;
    std::size_t start;
    bool is_last;
};

class TaskQueue {
public:
    TaskQueue(std::size_t capacity);
    ~TaskQueue();
    void push(Task task);
    bool pop(Task& out);
    void close();

private:
    std::mutex m_;
    std::condition_variable cv_not_empty_, cv_not_full_;
    std::deque<Task> q_;
    std::size_t cap_;
    bool closed_ = false;
};
