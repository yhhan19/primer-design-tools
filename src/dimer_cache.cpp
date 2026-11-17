#include "dimer_cache.hpp"
#include <random>

namespace {
std::string digest(const std::string& value) {
    std::uint64_t hash = 14695981039346656037ULL;
    for (unsigned char c : value) { hash ^= c; hash *= 1099511628211ULL; }
    std::ostringstream out;
    out << std::hex << hash;
    return out.str();
}
std::shared_future<double> ready(double value) {
    std::promise<double> p;
    p.set_value(value);
    return p.get_future().share();
}
}

DimerInteractionCache::DimerInteractionCache(std::string context,
                                             const std::string& directory)
    : context_(std::move(context)) {
    if (directory.empty()) return; // In-memory cache.
    directory_ = (fs::path(directory) / ("v1-" + digest(context_))).string();
    fs::create_directories(directory_);
    for (const auto& item : fs::directory_iterator(directory_)) {
        if (!item.is_regular_file() || item.path().extension() != ".tsv") continue;
        std::ifstream input(item.path(), std::ios::binary);
        std::string header;
        std::getline(input, header);
        std::ostringstream payload;
        payload << input.rdbuf();
        const auto data = payload.str();
        if (!input || header != "dpro-dimer-cache-v1\t" + digest(data))
            throw std::runtime_error("Invalid dimer cache shard: " + item.path().string());
        std::istringstream rows(data);
        std::string context_line;
        std::getline(rows, context_line);
        // Compare full context as well as its directory hash: hash collisions
        // can cause a miss, never use of another context's energies.
        if (context_line != context_) continue;
        std::string line;
        while (std::getline(rows, line)) {
            const auto tab = line.rfind('\t');
            if (tab == std::string::npos)
                throw std::runtime_error("Malformed dimer cache row: " + item.path().string());
            std::size_t consumed = 0;
            const auto token = line.substr(tab + 1);
            const double value = std::stod(token, &consumed);
            if (consumed != token.size() || !std::isfinite(value))
                throw std::runtime_error("Invalid dimer cache energy: " + item.path().string());
            const auto key = line.substr(0, tab);
            auto old = entries_.find(key);
            if (old != entries_.end() && old->second.get() != value)
                throw std::runtime_error("Conflicting dimer cache energies: " + item.path().string());
            entries_.emplace(key, ready(value));
        }
    }
}

double DimerInteractionCache::get(const std::string& first, const std::string& second,
    double min_c, double max_c, const std::function<double()>& calculate) {
    std::ostringstream key_stream;
    key_stream << std::setprecision(std::numeric_limits<double>::max_digits10)
               << min_c << '\t' << max_c << '\t' << first << '\t' << second;
    const auto key = key_stream.str();
    std::shared_future<double> future;
    std::shared_ptr<std::promise<double>> producer;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        auto found = entries_.find(key);
        if (found != entries_.end()) {
            ++stats_.reused;
            future = found->second;
        } else {
            producer = std::make_shared<std::promise<double>>();
            future = producer->get_future().share();
            entries_.emplace(key, future);
        }
    }
    if (!producer) return future.get();
    try {
        const double result = calculate();
        if (!std::isfinite(result)) throw std::runtime_error("Non-finite cached dimer energy");
        {
            std::lock_guard<std::mutex> lock(mutex_);
            ++stats_.computed;
            pending_[key] = result;
        }
        producer->set_value(result);
        return result;
    } catch (...) {
        producer->set_exception(std::current_exception());
        std::lock_guard<std::mutex> lock(mutex_);
        entries_.erase(key);
        throw;
    }
}

void DimerInteractionCache::flush() {
    std::lock_guard<std::mutex> writer(flush_mutex_);
    std::map<std::string, double> batch;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (directory_.empty()) { pending_.clear(); return; }
        batch.swap(pending_);
    }
    if (batch.empty()) return;
    fs::path staging;
    try {
        std::random_device random;
        do {
            staging = fs::path(directory_) / (".pending-" + std::to_string(random()) +
                "-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        } while (!fs::create_directory(staging));
        std::ostringstream payload;
        payload << context_ << '\n' << std::setprecision(std::numeric_limits<double>::max_digits10);
        for (const auto& entry : batch) payload << entry.first << '\t' << entry.second << '\n';
        const auto data = payload.str();
        const auto temporary = staging / "entries";
        std::ofstream output(temporary, std::ios::binary);
        output << "dpro-dimer-cache-v1\t" << digest(data) << '\n' << data;
        output.close();
        if (!output) throw std::runtime_error("Cannot write dimer cache: " + temporary.string());
        // Publish complete immutable shards only. Concurrent processes cannot
        // overwrite each other's entries; interrupted temporary files are ignored.
        fs::rename(temporary, fs::path(directory_) / (staging.filename().string() + ".tsv"));
        fs::remove(staging);
    } catch (...) {
        std::lock_guard<std::mutex> lock(mutex_);
        pending_.insert(batch.begin(), batch.end());
        throw;
    }
}

DimerInteractionCache::Stats DimerInteractionCache::stats() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return stats_;
}
