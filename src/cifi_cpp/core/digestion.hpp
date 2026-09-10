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
    // Minimum length of an *emitted* read, guaranteed for both R1 and R2.
    int min_frag_len = 60;
    // Drop the 5' site remnant from fragments that begin at a cut.
    bool strip_overhang = true;
    // Reverse complement R2 (after any strip). Off by default: fragments keep
    // the orientation they were sequenced in, as in the Pore-C convention.
    bool revcomp_r2 = false;
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
 * Extract the emitted span of each fragment.
 *
 * lead_trim is removed from fragments that begin at a cut site; the read's
 * leading fragment does not begin at one and so is never trimmed. Spans shorter
 * than min_emit_len are dropped, which makes min_emit_len a guarantee about the
 * emitted read rather than about the untrimmed source fragment.
 */
std::vector<std::pair<size_t, size_t>> extract_fragments(
    const std::string& sequence,
    const EnzymeInfo& enzyme,
    int min_emit_len,
    int lead_trim = 0
);

} // namespace cifi
