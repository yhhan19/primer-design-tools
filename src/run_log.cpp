#include "run_log.hpp"
#include <algorithm>
#include <chrono>
#include <csignal>
#include <cstdlib>
#include <deque>
#include <exception>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <vector>
#ifdef _WIN32
#include <io.h>
#include <windows.h>
#else
#include <sys/ioctl.h>
#include <unistd.h>
#endif

namespace {
RunLog* current_log = nullptr;
#ifndef _WIN32
volatile std::sig_atomic_t terminal_screen_active = 0;
void restore_terminal_on_signal(int signal) {
    if (terminal_screen_active) {
        constexpr char restore[] = "\033[?25h\033[?1049l";
        const auto ignored = ::write(STDOUT_FILENO, restore, sizeof(restore)-1);
        (void)ignored;
    }
    _exit(128+signal);
}
#endif
bool terminal_available() {
#ifdef _WIN32
    if (!_isatty(_fileno(stdout))) return false;
    const HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD mode = 0;
    return GetConsoleMode(handle, &mode) && SetConsoleMode(handle, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
#else
    const char* term = std::getenv("TERM");
    return isatty(STDOUT_FILENO) && (!term || std::string(term) != "dumb");
#endif
}
std::pair<std::size_t, std::size_t> terminal_size() {
#ifndef _WIN32
    winsize size{};
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &size) == 0 && size.ws_row && size.ws_col)
        return {std::max<std::size_t>(4, size.ws_row), std::max<std::size_t>(10, size.ws_col)};
#endif
    return {30, 100};
}
// Terminal text is clipped at a UTF-8 boundary, with control bytes removed.
std::string screen_text(const std::string& input, std::size_t columns) {
    std::string out;
    std::size_t used = 0;
    for (std::size_t i = 0; i < input.size() && used < columns;) {
        const auto c = static_cast<unsigned char>(input[i]);
        if (c < 32 || c == 127) { ++i; continue; }
        std::size_t bytes = c < 128 ? 1 : c < 224 ? 2 : c < 240 ? 3 : 4;
        bytes = std::min(bytes, input.size()-i);
        const std::size_t cells = bytes > 1 ? 2 : 1;
        if (used+cells > columns) break;
        out.append(input, i, bytes); i += bytes; used += cells;
    }
    return out;
}
}

struct RunLog::Impl {
    using Clock = std::chrono::steady_clock;
    struct Frame {
        std::size_t id, completed_start;
        std::string title, partial;
        std::deque<std::string> lines;
        Clock::time_point start = Clock::now();
    };
    struct Buffer : std::streambuf {
        Impl& owner;
        std::streambuf* original;
        Buffer(Impl& owner, std::streambuf* original) : owner(owner), original(original) {}
        std::streamsize xsputn(const char* p, std::streamsize n) override {
            owner.write(p, static_cast<std::size_t>(n), original); return n;
        }
        int_type overflow(int_type value) override {
            if (traits_type::eq_int_type(value, traits_type::eof())) return traits_type::not_eof(value);
            const char c = traits_type::to_char_type(value); owner.write(&c, 1, original); return value;
        }
        int sync() override {
            std::lock_guard<std::mutex> lock(owner.mutex);
            if (owner.file.is_open()) owner.file.flush();
            original->pubsync(); owner.render(false); return 0;
        }
    };
    std::streambuf *out_original = std::cout.rdbuf(), *err_original = std::cerr.rdbuf();
    Buffer out{*this, out_original}, err{*this, err_original};
    std::mutex mutex;
    std::ofstream file;
    std::string pending, path;
    bool live = false, alternate = false, finished = false;
#ifndef _WIN32
    using SignalHandler = void (*)(int);
    SignalHandler old_int = SIG_DFL, old_term = SIG_DFL;
    bool signals_installed = false;
#endif
    std::size_t next_id = 1;
    std::vector<Frame> frames;
    std::vector<std::string> completed;
    Clock::time_point last_render{};

    explicit Impl(bool allow_live) {
        frames.push_back({0, 0, "Run output", "", {}});
        live = allow_live && terminal_available();
        std::cout.rdbuf(&out); std::cerr.rdbuf(&err);
    }
    ~Impl() {
#ifndef _WIN32
        if (signals_installed) { std::signal(SIGINT, old_int); std::signal(SIGTERM, old_term); }
#endif
    }
    void terminal_write(const std::string& value) {
        out_original->sputn(value.data(), value.size()); out_original->pubsync();
    }
    void enter() {
        if (!alternate) {
#ifndef _WIN32
            if (!signals_installed) {
                old_int = std::signal(SIGINT, restore_terminal_on_signal);
                old_term = std::signal(SIGTERM, restore_terminal_on_signal);
                signals_installed = true;
            }
            terminal_screen_active = 1;
#endif
            terminal_write("\033[?1049h\033[?25l"); alternate = true;
        }
    }
    void leave() {
        if (alternate) {
            terminal_write("\033[?25h\033[?1049l"); alternate = false;
#ifndef _WIN32
            terminal_screen_active = 0;
#endif
        }
    }
    void write(const char* p, std::size_t n, std::streambuf* original) {
        std::lock_guard<std::mutex> lock(mutex);
        if (file.is_open()) file.write(p, static_cast<std::streamsize>(n));
        else pending.append(p, n);
        bool newline = false;
        auto& frame = frames.back();
        for (std::size_t i=0; i<n; ++i) {
            if (p[i] == '\n') {
                frame.lines.push_back(std::move(frame.partial)); frame.partial.clear();
                if (frame.lines.size() > 2000) frame.lines.pop_front();
                newline = true;
            } else if (p[i] == '\r') frame.partial.clear();
            else frame.partial += p[i];
        }
        if (!live || finished) original->sputn(p, static_cast<std::streamsize>(n));
        if (newline) { if (file.is_open()) file.flush(); render(false); }
    }
    void render(bool force) {
        if (!live || finished || !alternate) return;
        const auto now = Clock::now();
        if (!force && now-last_render < std::chrono::milliseconds(75)) return;
        last_render = now;
        const auto [height, width] = terminal_size();
        std::vector<std::string> display;
        display.push_back("dpro | full log: " + (path.empty() ? "initializing" : path));
        const std::size_t visible_completed = std::min(completed.size(), height>8 ? height/3 : 0);
        if (completed.size() > visible_completed && height>8)
            display.push_back(std::to_string(completed.size()-visible_completed) + " earlier blocks completed (in full log)");
        for (std::size_t i=completed.size()-visible_completed; i<completed.size(); ++i)
            display.push_back(completed[i]);
        std::string title = "Running: ";
        for (std::size_t i=1; i<frames.size(); ++i) {
            if (i>1) title += " / "; title += frames[i].title;
        }
        if (frames.size()==1) title += frames[0].title;
        display.push_back(title);
        display.push_back(std::string(std::min<std::size_t>(width-1, 100), '-'));
        const auto& frame = frames.back();
        const std::size_t available = height > display.size()+1 ? height-display.size()-1 : 0;
        const std::size_t count = frame.lines.size() + !frame.partial.empty();
        const std::size_t first = count > available ? count-available : 0;
        if (available) {
            for (std::size_t i=first; i<frame.lines.size(); ++i) display.push_back(frame.lines[i]);
            if (!frame.partial.empty()) display.push_back(frame.partial);
        }
        std::ostringstream screen;
        screen << "\033[H\033[2J";
        for (std::size_t i=0; i<display.size(); ++i) {
            if (i) screen << "\r\n";
            screen << screen_text(display[i], width-1);
        }
        terminal_write(screen.str());
    }
};

RunLog::RunLog(bool allow_live) : impl_(std::make_unique<Impl>(allow_live)) { current_log = this; }
RunLog::~RunLog() {
    try { finish(true); } catch (...) { impl_->leave(); }
    std::cout.rdbuf(impl_->out_original); std::cerr.rdbuf(impl_->err_original);
    current_log = nullptr;
}
RunLog* RunLog::active() { return current_log; }
void RunLog::configure(const std::string& path, const std::string& mode) {
    std::lock_guard<std::mutex> lock(impl_->mutex);
    if (impl_->file.is_open()) throw std::runtime_error("Run log is already open");
    impl_->file.open(path, std::ios::out | std::ios::trunc);
    if (!impl_->file) throw std::runtime_error("Cannot write full run log: " + path);
    impl_->path = path;
    impl_->file << impl_->pending; impl_->file.flush();
    const bool live = mode == "live" || (mode == "auto" && terminal_available());
    if (impl_->live && !live) { impl_->leave(); impl_->terminal_write(impl_->pending); }
    impl_->pending.clear(); impl_->live = live;
    if (live) { impl_->enter(); impl_->render(true); }
}
std::size_t RunLog::begin(const std::string& title) {
    std::lock_guard<std::mutex> lock(impl_->mutex);
    impl_->frames.back().lines.clear(); impl_->frames.back().partial.clear();
    const auto id = impl_->next_id++;
    impl_->frames.push_back({id, impl_->completed.size(), title, "", {}});
    impl_->render(true); return id;
}
void RunLog::end(std::size_t id, const std::string& summary, bool failed) {
    std::lock_guard<std::mutex> lock(impl_->mutex);
    if (impl_->frames.size()<2 || impl_->frames.back().id != id) return;
    const auto frame = std::move(impl_->frames.back()); impl_->frames.pop_back();
    impl_->completed.resize(frame.completed_start);
    const double seconds = std::chrono::duration<double>(Impl::Clock::now()-frame.start).count();
    std::ostringstream line;
    line << (failed ? "[FAILED] " : "[done] ") << frame.title;
    if (!summary.empty()) line << ": " << summary;
    line << " (" << std::fixed << std::setprecision(1) << seconds << " s)";
    impl_->completed.push_back(line.str());
    if (impl_->file.is_open()) { impl_->file << line.str() << '\n'; impl_->file.flush(); }
    else impl_->pending += line.str()+"\n";
    if (!impl_->live) impl_->terminal_write(line.str()+"\n");
    impl_->render(true);
}
void RunLog::archive_error_if_needed() {
    if (!impl_->file.is_open()) {
        try { configure("dpro-error.log", "auto"); } catch (...) {}
    }
}
void RunLog::finish(bool success, const std::string& message) {
    std::lock_guard<std::mutex> lock(impl_->mutex);
    if (impl_->finished) return;
    if (impl_->live) {
        impl_->leave();
        for (const auto& line : impl_->completed) impl_->terminal_write(line+"\n");
        // Keep the final report visible after restoring the user's terminal.
        const auto& frame = impl_->frames.back();
        const std::size_t start = frame.lines.size()>25 ? frame.lines.size()-25 : 0;
        for (std::size_t i=start; i<frame.lines.size(); ++i) impl_->terminal_write(frame.lines[i]+"\n");
        if (!frame.partial.empty()) impl_->terminal_write(frame.partial+"\n");
    }
    if (!message.empty()) {
        const std::string line = (success ? "" : "Run failed: ")+message+"\n";
        if (impl_->file.is_open()) impl_->file << line;
        else impl_->pending += line;
        impl_->terminal_write(line);
    }
    if (impl_->file.is_open()) {
        const std::string line = "Full run log: " + impl_->path + "\n";
        impl_->file << line; impl_->file.flush(); impl_->terminal_write(line);
    }
    impl_->finished = true;
}
RunLogBlock::RunLogBlock(const std::string& title)
    : log_(RunLog::active()), id_(log_ ? log_->begin(title) : 0), exceptions_(std::uncaught_exceptions()) {}
RunLogBlock::~RunLogBlock() {
    if (log_) log_->end(id_, "", std::uncaught_exceptions()>exceptions_);
}
void RunLogBlock::complete(const std::string& summary) {
    if (log_) { log_->end(id_, summary, false); log_=nullptr; }
}
void configure_run_log(const std::string& path, const std::string& mode) {
    if (auto* log = RunLog::active()) log->configure(path, mode);
}
