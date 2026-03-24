#include "task.hpp"

TaskQueue::TaskQueue(std::size_t capacity) : cap_(capacity) {}

TaskQueue::~TaskQueue() {}

void TaskQueue::push(Task task) {
    std::unique_lock<std::mutex> lk(m_);
    cv_not_full_.wait(lk, [&]{ return closed_ || q_.size() < cap_; });
    if (closed_) return;
    q_.push_back(std::move(task));
    cv_not_empty_.notify_one();
}

bool TaskQueue::pop(Task& out) {
    std::unique_lock<std::mutex> lk(m_);
    cv_not_empty_.wait(lk, [&]{ return closed_ || !q_.empty(); });
    if (q_.empty()) return false;
    out = std::move(q_.front());
    q_.pop_front();
    cv_not_full_.notify_one();
    return true;
}

void TaskQueue::close() {
    std::lock_guard<std::mutex> lk(m_);
    closed_ = true;
    cv_not_empty_.notify_all();
    cv_not_full_.notify_all();
}
