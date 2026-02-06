#include "statistics.hpp"
#include <algorithm>
#include <cmath>

namespace cifi {

Statistics::Statistics(bool fast_mode, int bin_size)
    : fast_mode_(fast_mode), bin_size_(bin_size) {}

void Statistics::add(int value) {
    count_++;
    sum_ += value;
    if (value < min_) min_ = value;
    if (value > max_) max_ = value;

    if (fast_mode_) {
        histogram_[value / bin_size_]++;
    } else {
        values_.push_back(value);
    }
}

void Statistics::clear() {
    count_ = 0;
    sum_ = 0;
    min_ = INT32_MAX;
    max_ = INT32_MIN;
    values_.clear();
    histogram_.clear();
}

double Statistics::mean() const {
    return count_ > 0 ? static_cast<double>(sum_) / count_ : 0.0;
}

double Statistics::median() const {
    return percentile(0.5);
}

double Statistics::percentile(double p) const {
    if (count_ == 0) return 0.0;

    // Normalize: accept both 0-1 and 0-100 ranges
    if (p > 1.0) p /= 100.0;
    if (p < 0.0) p = 0.0;
    if (p > 1.0) p = 1.0;

    if (fast_mode_) {
        return percentile_from_histogram(p);
    }

    // Exact mode: sort and index
    std::vector<int> sorted = values_;
    std::sort(sorted.begin(), sorted.end());

    size_t idx = static_cast<size_t>(p * (sorted.size() - 1));
    return sorted[idx];
}

double Statistics::percentile_from_histogram(double p) const {
    if (histogram_.empty()) return 0.0;

    // Get sorted bins
    std::vector<std::pair<int, uint64_t>> bins(histogram_.begin(), histogram_.end());
    std::sort(bins.begin(), bins.end());

    // Find target count
    uint64_t target = static_cast<uint64_t>(p * count_);
    uint64_t cumulative = 0;

    for (const auto& [bin, cnt] : bins) {
        cumulative += cnt;
        if (cumulative >= target) {
            // Interpolate within bin
            double bin_start = bin * bin_size_;
            double bin_end = bin_start + bin_size_;
            return (bin_start + bin_end) / 2.0;
        }
    }

    return bins.back().first * bin_size_ + bin_size_ / 2.0;
}

std::vector<std::pair<int, uint64_t>> Statistics::get_histogram() const {
    if (fast_mode_) {
        std::vector<std::pair<int, uint64_t>> result(histogram_.begin(), histogram_.end());
        std::sort(result.begin(), result.end());
        return result;
    }

    // Build histogram from raw values
    phmap::flat_hash_map<int, uint64_t> hist;
    for (int v : values_) {
        hist[v / bin_size_]++;
    }

    std::vector<std::pair<int, uint64_t>> result(hist.begin(), hist.end());
    std::sort(result.begin(), result.end());
    return result;
}

} // namespace cifi
