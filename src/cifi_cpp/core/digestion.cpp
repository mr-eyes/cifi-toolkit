#include "digestion.hpp"

namespace cifi {

std::vector<std::pair<size_t, size_t>> extract_fragments(
    const std::string& sequence,
    const EnzymeInfo& enzyme,
    int min_emit_len,
    int lead_trim
) {
    auto sites = find_all_degenerate(sequence, enzyme.site);

    // Build cut positions
    std::vector<size_t> cuts;
    cuts.push_back(0);
    for (size_t pos : sites) {
        cuts.push_back(pos + enzyme.cut_offset);
    }
    cuts.push_back(sequence.length());

    // Take each fragment's emitted span, then filter on that length. Only a
    // fragment that begins at a cut carries the site remnant; the read's
    // leading fragment starts at position 0 and so keeps its full length.
    std::vector<std::pair<size_t, size_t>> fragments;
    for (size_t i = 0; i < cuts.size() - 1; i++) {
        size_t start = cuts[i];
        size_t end = cuts[i + 1];
        if (start > 0) {
            start += static_cast<size_t>(lead_trim);
        }
        if (end > start && static_cast<int>(end - start) >= min_emit_len) {
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

    // The site remnant belongs to the fragment, not to whichever slot the
    // fragment lands in, so trim once here and let both mates read the same
    // spans. min_frag_len then bounds the emitted read directly.
    int lead_trim = config.strip_overhang ? config.enzyme.overhang_length() : 0;

    auto fragments = extract_fragments(sequence, config.enzyme,
                                       config.min_frag_len, lead_trim);

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

    // Generate ALL pairs (n choose 2)
    for (size_t i = 0; i < fragments.size(); i++) {
        for (size_t j = i + 1; j < fragments.size(); j++) {
            const auto& f1 = fragments[i];
            const auto& f2 = fragments[j];

            std::string seq1 = sequence.substr(f1.first, f1.second - f1.first);
            std::string qual1 = quality.substr(f1.first, f1.second - f1.first);
            std::string seq2 = sequence.substr(f2.first, f2.second - f2.first);
            std::string qual2 = quality.substr(f2.first, f2.second - f2.first);

            // Fragments arrive already trimmed, so R2 differs from R1 only by
            // the optional reverse complement.
            std::string r2_seq, r2_qual;
            if (config.revcomp_r2) {
                r2_seq = revcomp(seq2);
                r2_qual.assign(qual2.rbegin(), qual2.rend());
            } else {
                r2_seq = std::move(seq2);
                r2_qual = std::move(qual2);
            }

            // Both mates carry one name: R1/R2 FASTQ has no flags, so the name
            // is the only thing that identifies a pair. A "/1" or "/2" suffix
            // would make the two files disagree on every pair.
            std::string pair_name = name + "_" + std::to_string(i) + "_" + std::to_string(j - i - 1);

            out_r1.write(pair_name, seq1, qual1);
            out_r2.write(pair_name, r2_seq, r2_qual);

            result.pairs_written++;
        }
    }

    return true;
}

} // namespace cifi
