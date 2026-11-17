#pragma once
#include <cstddef>
#include <memory>
#include <string>

// Complete plain-text archive plus a live, collapsible terminal view.
class RunLog {
public:
    explicit RunLog(bool allow_live = true);
    ~RunLog();
    RunLog(const RunLog&) = delete;
    RunLog& operator=(const RunLog&) = delete;
    void configure(const std::string& path, const std::string& console_mode);
    void finish(bool success, const std::string& message = "");
    void archive_error_if_needed();
    static RunLog* active();
    std::size_t begin(const std::string& title);
    void end(std::size_t id, const std::string& summary, bool failed);
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class RunLogBlock {
public:
    explicit RunLogBlock(const std::string& title);
    ~RunLogBlock();
    RunLogBlock(const RunLogBlock&) = delete;
    RunLogBlock& operator=(const RunLogBlock&) = delete;
    void complete(const std::string& summary = "");
private:
    RunLog* log_;
    std::size_t id_;
    int exceptions_;
};

void configure_run_log(const std::string& path, const std::string& console_mode);
