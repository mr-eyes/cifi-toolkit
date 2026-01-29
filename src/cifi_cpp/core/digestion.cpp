#include "digestion.hpp"

namespace cifi {

std::vector<std::pair<size_t, size_t>> extract_fragments(
    const std::string& sequence,
    const EnzymeInfo& enzyme,
    int min_frag_len
) {
    auto sites = find_all_degenerate(sequence, enzyme.site);

    // Build cut positions
    std::vector<size_t> cuts;
    cuts.push_back(0);
    for (size_t pos : sites) {
        cuts.push_back(pos + enzyme.cut_offset);
    }
    cuts.push_back(sequence.length());

    // Extract fragments meeting length requirement
    std::vector<std::pair<size_t, size_t>> fragments;
    for (size_t i = 0; i < cuts.size() - 1; i++) {
        size_t start = cuts[i];
        size_t end = cuts[i + 1];
        if (end > start && static_cast<int>(end - start) >= min_frag_len) {
            fragments.push_back({start, end});
        }
    }

    return fragments;
}

bool process_single_read(
    const std::string& name,
    const std::string& sequence,
    const std::string& quality,
    const ProcessingConfig& config,
    FastqWriter& out_r1,
    FastqWriter& out_r2,
    ProcessingResult& result
) {
    // Find sites (supports degenerate IUPAC bases)
    auto sites = find_all_degenerate(sequence, config.enzyme.site);

    // Early exit: not enough sites
    if (static_cast<int>(sites.size()) < config.min_fragments - 1) {
        result.reads_skipped++;
        result.filtered_few_sites++;
        return false;
    }

    // Extract fragments
    auto fragments = extract_fragments(sequence, config.enzyme, config.min_frag_len);

    // Check fragment count after length filtering
    if (static_cast<int>(fragments.size()) < config.min_fragments) {
        result.reads_skipped++;
        result.filtered_short_frags++;
        return false;
    }

    // Record stats for passing reads only
    result.sites_per_read_stats.add(static_cast<int>(sites.size()));
    for (const auto& [start, end] : fragments) {
        result.frag_length_stats.add(static_cast<int>(end - start));
    }

    result.reads_out++;
    result.total_frags += fragments.size();

    // Calculate overhang for R2 stripping
    int overhang = config.enzyme.overhang_length();

    // Generate ALL pairs (n choose 2)
    for (size_t i = 0; i < fragments.size(); i++) {
        for (size_t j = i + 1; j < fragments.size(); j++) {
            const auto& f1 = fragments[i];
            const auto& f2 = fragments[j];

            std::string seq1 = sequence.substr(f1.first, f1.second - f1.first);
            std::string qual1 = quality.substr(f1.first, f1.second - f1.first);
            std::string seq2 = sequence.substr(f2.first, f2.second - f2.first);
            std::string qual2 = quality.substr(f2.first, f2.second - f2.first);

            std::string r2_seq, r2_qual;
            if (config.strip_overhang) {
                if (static_cast<int>(seq2.length()) > overhang) {
                    r2_seq = seq2.substr(overhang);
                    r2_qual = qual2.substr(overhang);
                } else {
                    r2_seq = seq2;
                    r2_qual = qual2;
                }
            } else {
                r2_seq = revcomp(seq2);
                r2_qual = std::string(qual2.rbegin(), qual2.rend());
            }

            // Build read names
            std::string r1_name = name + "_" + std::to_string(i) + "_" + std::to_string(j - i - 1) + "/1";
            std::string r2_name = name + "_" + std::to_string(i) + "_" + std::to_string(j - i - 1) + "/2";

            out_r1.write(r1_name, seq1, qual1);
            out_r2.write(r2_name, r2_seq, r2_qual);

            result.pairs_written++;
        }
    }

    return true;
}

} // namespace cifi
