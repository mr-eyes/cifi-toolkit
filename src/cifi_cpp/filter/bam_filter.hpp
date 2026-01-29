// successor to https://github.com/mydennislab/2024-sep-mapqfilter
#pragma once

#include <string>
#include <cstdint>

namespace cifi {

// Result struct for Python bindings
struct FilterResult {
    uint64_t total_reads = 0;
    uint64_t total_pairs = 0;
    uint64_t passed_pairs = 0;
    uint64_t failed_unpaired = 0;
    uint64_t failed_unmapped = 0;
    uint64_t failed_mapq = 0;
    uint64_t failed_mate_not_found = 0;
};

// Main filter function
FilterResult filter_bam(
    const std::string& input_path,
    const std::string& output_path,
    int mapq_threshold,
    int threads = 4
);

}  // namespace cifi
