#pragma once

#include "enzyme.hpp"
#include "../stats/statistics.hpp"
#include "../io/writer.hpp"
#include <string>
#include <vector>
#include <memory>

namespace cifi {

struct ProcessingConfig {
    EnzymeInfo enzyme;
    int min_fragments = 3;
    int min_frag_len = 20;
    bool strip_overhang = true;
    bool fast_mode = false;
};

struct ProcessingResult {
    uint64_t reads_in = 0;
    uint64_t reads_out = 0;
    uint64_t reads_skipped = 0;
    uint64_t pairs_written = 0;
    uint64_t total_frags = 0;

    // Filtering reason counters
    uint64_t filtered_few_sites = 0;    // Reads with < min_fragments sites
    uint64_t filtered_short_frags = 0;  // Reads where all fragments were too short after length filtering

    Statistics frag_length_stats;
    Statistics sites_per_read_stats;

    ProcessingResult(bool fast_mode = false)
        : frag_length_stats(fast_mode, 100)
        , sites_per_read_stats(fast_mode, 1) {}
};

/**
 * Process a single read: digest and write all pairwise contacts.
 * Returns true if read passed filters and was processed.
 */
bool process_single_read(
    const std::string& name,
    const std::string& sequence,
    const std::string& quality,
    const ProcessingConfig& config,
    FastqWriter& out_r1,
    FastqWriter& out_r2,
    ProcessingResult& result
);

/**
 * Extract fragments from a sequence given enzyme info.
 */
std::vector<std::pair<size_t, size_t>> extract_fragments(
    const std::string& sequence,
    const EnzymeInfo& enzyme,
    int min_frag_len
);

} // namespace cifi
