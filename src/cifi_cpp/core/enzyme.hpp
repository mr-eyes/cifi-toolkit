#pragma once

#include <string>
#include <vector>
#include <utility>
#include <optional>

namespace cifi {

struct EnzymeInfo {
    std::string name;
    std::string site;
    int cut_offset;

    int overhang_length() const;
};

/**
 * Get enzyme info by name.
 * Returns nullopt if enzyme not found.
 */
std::optional<EnzymeInfo> get_enzyme(const std::string& name);

/**
 * Get list of all available enzyme names.
 */
std::vector<std::string> list_enzymes();

/**
 * Find all occurrences of pattern in text.
 */
std::vector<size_t> find_all(const std::string& text, const std::string& pattern);

/**
 * Find all occurrences of pattern in text, supporting IUPAC degenerate bases.
 * Supported codes: N (any), R (A/G), Y (C/T), W (A/T), S (C/G),
 *                  M (A/C), K (G/T), B (C/G/T), D (A/G/T), H (A/C/T), V (A/C/G)
 */
std::vector<size_t> find_all_degenerate(const std::string& text, const std::string& pattern);

/**
 * Check if a pattern contains degenerate IUPAC bases.
 */
bool has_degenerate_bases(const std::string& pattern);

/**
 * Reverse complement a DNA sequence.
 */
std::string revcomp(const std::string& seq);

/**
 * Reverse complement a degenerate DNA sequence (handles IUPAC codes).
 */
std::string revcomp_degenerate(const std::string& seq);

} // namespace cifi
