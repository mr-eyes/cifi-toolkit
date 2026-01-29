// cifi - Toolkit for downstream processing of CiFi long reads.
// https://dennislab.org/cifi

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>

#include <zlib.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>

#include "core/enzyme.hpp"
#include "core/digestion.hpp"
#include "stats/statistics.hpp"
#include "io/writer.hpp"
#include "filter/bam_filter.hpp"

extern "C" {
#include "kseq.h"
}

KSEQ_INIT(gzFile, gzread)

namespace nb = nanobind;

// ============================================================================
// Helper functions
// ============================================================================

bool is_bam_file(const std::string& path) {
    size_t dot = path.rfind('.');
    if (dot != std::string::npos) {
        std::string ext = path.substr(dot);
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        return (ext == ".bam" || ext == ".sam" || ext == ".cram");
    }
    return false;
}

// ============================================================================
// QC: Single-enzyme analysis
// ============================================================================

struct SingleEnzymeQCResult {
    // Read statistics
    uint64_t reads_analyzed = 0;
    uint64_t total_bases = 0;
    double avg_read_length = 0;
    double median_read_length = 0;
    uint64_t min_read_length = 0;
    uint64_t max_read_length = 0;
    double gc_content = 0;

    // Site statistics
    uint64_t total_sites = 0;
    double sites_per_read_mean = 0;
    double sites_per_read_median = 0;
    int sites_per_read_min = 0;
    int sites_per_read_max = 0;

    // Yield estimates
    uint64_t reads_passing = 0;
    double pass_rate = 0;
    uint64_t est_total_fragments = 0;
    uint64_t est_total_pairs = 0;
    double avg_fragments_per_read = 0;
    double avg_pairs_per_read = 0;

    // Fragment size statistics
    double frag_size_mean = 0;
    double frag_size_median = 0;
    int frag_size_min = 0;
    int frag_size_max = 0;

    // Histograms (bin edges and counts)
    std::vector<double> read_length_hist_bins;
    std::vector<uint64_t> read_length_hist_counts;
    std::vector<double> sites_hist_bins;
    std::vector<uint64_t> sites_hist_counts;
    std::vector<double> frag_size_hist_bins;
    std::vector<uint64_t> frag_size_hist_counts;
};

static double compute_median(std::vector<int>& values) {
    if (values.empty()) return 0;
    std::sort(values.begin(), values.end());
    size_t n = values.size();
    if (n % 2 == 0) {
        return (values[n/2 - 1] + values[n/2]) / 2.0;
    }
    return values[n/2];
}

static void make_histogram(const std::vector<int>& data, int num_bins,
                           std::vector<double>& bins, std::vector<uint64_t>& counts) {
    bins.clear();
    counts.clear();

    if (data.empty()) return;

    int min_val = *std::min_element(data.begin(), data.end());
    int max_val = *std::max_element(data.begin(), data.end());

    if (min_val == max_val) {
        bins.push_back(min_val);
        counts.push_back(data.size());
        return;
    }

    double bin_width = static_cast<double>(max_val - min_val) / num_bins;
    bins.resize(num_bins + 1);
    counts.resize(num_bins, 0);

    for (int i = 0; i <= num_bins; ++i) {
        bins[i] = min_val + i * bin_width;
    }

    for (int val : data) {
        int idx = std::min(static_cast<int>((val - min_val) / bin_width), num_bins - 1);
        counts[idx]++;
    }
}

SingleEnzymeQCResult run_qc_analysis_custom(
    const std::string& input_path,
    const std::string& site,
    int cut_offset,
    int num_reads,
    int min_sites = 2
) {
    SingleEnzymeQCResult result;

    std::vector<int> read_lengths;
    std::vector<int> sites_per_read;
    std::vector<int> fragment_sizes;
    uint64_t gc_count = 0;
    uint64_t total_bases = 0;
    uint64_t min_len = UINT64_MAX;
    uint64_t max_len = 0;

    auto process_sequence = [&](const std::string& sequence) {
        size_t seq_len = sequence.length();
        result.reads_analyzed++;
        total_bases += seq_len;
        read_lengths.push_back(static_cast<int>(seq_len));

        if (seq_len < min_len) min_len = seq_len;
        if (seq_len > max_len) max_len = seq_len;

        for (char c : sequence) {
            if (c == 'G' || c == 'C' || c == 'g' || c == 'c') gc_count++;
        }

        auto positions = cifi::find_all_degenerate(sequence, site);
        int n_sites = static_cast<int>(positions.size());
        sites_per_read.push_back(n_sites);
        result.total_sites += n_sites;

        if (!positions.empty()) {
            if (positions[0] > 0) {
                int frag_len = static_cast<int>(positions[0]) + cut_offset;
                if (frag_len > 0) fragment_sizes.push_back(frag_len);
            }
            for (size_t i = 0; i + 1 < positions.size(); ++i) {
                int frag_len = static_cast<int>(positions[i + 1] - positions[i]);
                if (frag_len > 0) fragment_sizes.push_back(frag_len);
            }
            size_t last_end = positions.back() + cut_offset;
            if (last_end < seq_len) {
                int frag_len = static_cast<int>(seq_len - last_end);
                if (frag_len > 0) fragment_sizes.push_back(frag_len);
            }
        } else {
            fragment_sizes.push_back(static_cast<int>(seq_len));
        }
    };

    // Process input file
    if (is_bam_file(input_path)) {
        htsFile *fp = hts_open(input_path.c_str(), "r");
        if (!fp) throw std::runtime_error("Cannot open: " + input_path);

        sam_hdr_t *hdr = sam_hdr_read(fp);
        if (!hdr) {
            hts_close(fp);
            throw std::runtime_error("Cannot read header: " + input_path);
        }

        bam1_t *b = bam_init1();

        while (sam_read1(fp, hdr, b) >= 0 && (num_reads == 0 || result.reads_analyzed < static_cast<uint64_t>(num_reads))) {
            if (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;

            std::string sequence;
            uint8_t *seq_ptr = bam_get_seq(b);
            sequence.reserve(b->core.l_qseq);
            for (int i = 0; i < b->core.l_qseq; i++) {
                sequence += seq_nt16_str[bam_seqi(seq_ptr, i)];
            }

            if (b->core.flag & BAM_FREVERSE) {
                sequence = cifi::revcomp(sequence);
            }

            process_sequence(sequence);
        }

        bam_destroy1(b);
        sam_hdr_destroy(hdr);
        hts_close(fp);
    } else {
        gzFile fp = gzopen(input_path.c_str(), "r");
        if (!fp) throw std::runtime_error("Cannot open: " + input_path);
        kseq_t *seq = kseq_init(fp);

        while (kseq_read(seq) >= 0 && (num_reads == 0 || result.reads_analyzed < static_cast<uint64_t>(num_reads))) {
            std::string sequence(seq->seq.s, seq->seq.l);
            process_sequence(sequence);
        }

        kseq_destroy(seq);
        gzclose(fp);
    }

    // Compute statistics
    if (result.reads_analyzed == 0) return result;

    result.total_bases = total_bases;
    result.min_read_length = (min_len == UINT64_MAX) ? 0 : min_len;
    result.max_read_length = max_len;
    result.avg_read_length = static_cast<double>(total_bases) / result.reads_analyzed;
    result.gc_content = total_bases > 0 ? 100.0 * gc_count / total_bases : 0;
    result.median_read_length = compute_median(read_lengths);

    std::vector<int> sorted_sites = sites_per_read;
    std::sort(sorted_sites.begin(), sorted_sites.end());
    size_t n = sorted_sites.size();

    result.sites_per_read_mean = static_cast<double>(result.total_sites) / n;
    result.sites_per_read_median = sorted_sites[n / 2];
    result.sites_per_read_min = sorted_sites.front();
    result.sites_per_read_max = sorted_sites.back();

    for (int s : sites_per_read) {
        if (s >= min_sites) {
            result.reads_passing++;
            int n_frags = s + 1;
            result.est_total_fragments += n_frags;
            result.est_total_pairs += (n_frags * (n_frags - 1)) / 2;
        }
    }
    result.pass_rate = 100.0 * result.reads_passing / result.reads_analyzed;
    result.avg_fragments_per_read = result.reads_passing > 0 ?
        static_cast<double>(result.est_total_fragments) / result.reads_passing : 0;
    result.avg_pairs_per_read = result.reads_passing > 0 ?
        static_cast<double>(result.est_total_pairs) / result.reads_passing : 0;

    if (!fragment_sizes.empty()) {
        std::vector<int> sorted_frags = fragment_sizes;
        std::sort(sorted_frags.begin(), sorted_frags.end());
        size_t nf = sorted_frags.size();

        double frag_sum = 0;
        for (int fs : fragment_sizes) frag_sum += fs;

        result.frag_size_mean = frag_sum / nf;
        result.frag_size_median = sorted_frags[nf / 2];
        result.frag_size_min = sorted_frags.front();
        result.frag_size_max = sorted_frags.back();

        make_histogram(fragment_sizes, 50, result.frag_size_hist_bins, result.frag_size_hist_counts);
    }

    make_histogram(read_lengths, 50, result.read_length_hist_bins, result.read_length_hist_counts);
    int sites_bins = std::min(result.sites_per_read_max + 1, 50);
    make_histogram(sites_per_read, sites_bins, result.sites_hist_bins, result.sites_hist_counts);

    return result;
}

// ============================================================================
// Read Processing
// ============================================================================

// Helper: read BAM sequence into a read for processing
static void process_bam_reads(
    const std::string& input_path,
    const cifi::ProcessingConfig& config,
    cifi::FastqWriter& writer_r1,
    cifi::FastqWriter& writer_r2,
    cifi::ProcessingResult& result
) {
    htsFile *fp = hts_open(input_path.c_str(), "r");
    if (!fp) throw std::runtime_error("Cannot open: " + input_path);

    sam_hdr_t *hdr = sam_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        throw std::runtime_error("Cannot read header: " + input_path);
    }

    bam1_t *b = bam_init1();

    while (sam_read1(fp, hdr, b) >= 0) {
        result.reads_in++;

        if (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
            continue;
        }

        std::string name(bam_get_qname(b));

        std::string sequence;
        uint8_t *seq_ptr = bam_get_seq(b);
        sequence.reserve(b->core.l_qseq);
        for (int i = 0; i < b->core.l_qseq; i++) {
            sequence += seq_nt16_str[bam_seqi(seq_ptr, i)];
        }

        std::string quality;
        uint8_t *qual_ptr = bam_get_qual(b);
        if (qual_ptr[0] != 0xff) {
            quality.reserve(b->core.l_qseq);
            for (int i = 0; i < b->core.l_qseq; i++) {
                quality += (char)(qual_ptr[i] + 33);
            }
        } else {
            quality = std::string(b->core.l_qseq, 'I');
        }

        if (b->core.flag & BAM_FREVERSE) {
            sequence = cifi::revcomp(sequence);
            std::reverse(quality.begin(), quality.end());
        }

        cifi::process_single_read(name, sequence, quality, config,
                                  writer_r1, writer_r2, result);
    }

    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    hts_close(fp);
}

// Helper: read FASTQ sequence into a read for processing
static void process_fastq_reads(
    const std::string& input_path,
    const cifi::ProcessingConfig& config,
    cifi::FastqWriter& writer_r1,
    cifi::FastqWriter& writer_r2,
    cifi::ProcessingResult& result
) {
    gzFile fp = gzopen(input_path.c_str(), "r");
    if (!fp) throw std::runtime_error("Cannot open: " + input_path);
    kseq_t *seq = kseq_init(fp);

    while (kseq_read(seq) >= 0) {
        result.reads_in++;
        std::string name(seq->name.s, seq->name.l);
        std::string sequence(seq->seq.s, seq->seq.l);
        std::string quality(seq->qual.l ? std::string(seq->qual.s, seq->qual.l)
                                        : std::string(seq->seq.l, 'I'));

        cifi::process_single_read(name, sequence, quality, config,
                                  writer_r1, writer_r2, result);
    }

    kseq_destroy(seq);
    gzclose(fp);
}

// Process reads with a named enzyme
cifi::ProcessingResult process_reads(
    const std::string& input_path,
    const std::string& output_r1,
    const std::string& output_r2,
    const std::string& enzyme_name,
    int min_fragments = 3,
    int min_frag_len = 20,
    bool strip_overhang = true,
    bool gzip_output = false,
    bool fast_mode = false
) {
    auto enzyme_opt = cifi::get_enzyme(enzyme_name);
    if (!enzyme_opt) {
        throw std::runtime_error("Unknown enzyme: " + enzyme_name);
    }

    cifi::ProcessingConfig config;
    config.enzyme = *enzyme_opt;
    config.min_fragments = min_fragments;
    config.min_frag_len = min_frag_len;
    config.strip_overhang = strip_overhang;
    config.fast_mode = fast_mode;

    cifi::ProcessingResult result(fast_mode);
    auto writer_r1 = cifi::make_writer(output_r1, gzip_output);
    auto writer_r2 = cifi::make_writer(output_r2, gzip_output);

    if (is_bam_file(input_path)) {
        process_bam_reads(input_path, config, *writer_r1, *writer_r2, result);
    } else {
        process_fastq_reads(input_path, config, *writer_r1, *writer_r2, result);
    }

    writer_r1->close();
    writer_r2->close();
    return result;
}

// Process reads with a custom enzyme site
cifi::ProcessingResult process_reads_custom(
    const std::string& input_path,
    const std::string& output_r1,
    const std::string& output_r2,
    const std::string& site,
    int cut_offset,
    int min_fragments = 3,
    int min_frag_len = 20,
    bool strip_overhang = true,
    bool gzip_output = false,
    bool fast_mode = false
) {
    cifi::EnzymeInfo enzyme{"Custom", site, cut_offset};

    cifi::ProcessingConfig config;
    config.enzyme = enzyme;
    config.min_fragments = min_fragments;
    config.min_frag_len = min_frag_len;
    config.strip_overhang = strip_overhang;
    config.fast_mode = fast_mode;

    cifi::ProcessingResult result(fast_mode);
    auto writer_r1 = cifi::make_writer(output_r1, gzip_output);
    auto writer_r2 = cifi::make_writer(output_r2, gzip_output);

    if (is_bam_file(input_path)) {
        process_bam_reads(input_path, config, *writer_r1, *writer_r2, result);
    } else {
        process_fastq_reads(input_path, config, *writer_r1, *writer_r2, result);
    }

    writer_r1->close();
    writer_r2->close();
    return result;
}

// ============================================================================
// nanobind module definition
// ============================================================================

NB_MODULE(_core, m) {
    m.doc() = "cifi - Accurate long-read chromosome conformation capture";

    // Statistics class
    nb::class_<cifi::Statistics>(m, "Statistics")
        .def("count", &cifi::Statistics::count)
        .def("sum", &cifi::Statistics::sum)
        .def("mean", &cifi::Statistics::mean)
        .def("median", &cifi::Statistics::median)
        .def("min", &cifi::Statistics::min)
        .def("max", &cifi::Statistics::max)
        .def("percentile", &cifi::Statistics::percentile)
        .def("is_fast_mode", &cifi::Statistics::is_fast_mode)
        .def("values", &cifi::Statistics::values);

    // ProcessingResult
    nb::class_<cifi::ProcessingResult>(m, "ProcessingResult")
        .def_ro("reads_in", &cifi::ProcessingResult::reads_in)
        .def_ro("reads_out", &cifi::ProcessingResult::reads_out)
        .def_ro("pairs_written", &cifi::ProcessingResult::pairs_written)
        .def_ro("total_frags", &cifi::ProcessingResult::total_frags)
        .def_ro("reads_skipped", &cifi::ProcessingResult::reads_skipped)
        .def_ro("filtered_few_sites", &cifi::ProcessingResult::filtered_few_sites)
        .def_ro("filtered_short_frags", &cifi::ProcessingResult::filtered_short_frags)
        .def_ro("frag_length_stats", &cifi::ProcessingResult::frag_length_stats)
        .def_ro("sites_per_read_stats", &cifi::ProcessingResult::sites_per_read_stats);

    // SingleEnzymeQCResult
    nb::class_<SingleEnzymeQCResult>(m, "SingleEnzymeQCResult")
        .def_ro("reads_analyzed", &SingleEnzymeQCResult::reads_analyzed)
        .def_ro("total_bases", &SingleEnzymeQCResult::total_bases)
        .def_ro("avg_read_length", &SingleEnzymeQCResult::avg_read_length)
        .def_ro("median_read_length", &SingleEnzymeQCResult::median_read_length)
        .def_ro("min_read_length", &SingleEnzymeQCResult::min_read_length)
        .def_ro("max_read_length", &SingleEnzymeQCResult::max_read_length)
        .def_ro("gc_content", &SingleEnzymeQCResult::gc_content)
        .def_ro("total_sites", &SingleEnzymeQCResult::total_sites)
        .def_ro("sites_per_read_mean", &SingleEnzymeQCResult::sites_per_read_mean)
        .def_ro("sites_per_read_median", &SingleEnzymeQCResult::sites_per_read_median)
        .def_ro("sites_per_read_min", &SingleEnzymeQCResult::sites_per_read_min)
        .def_ro("sites_per_read_max", &SingleEnzymeQCResult::sites_per_read_max)
        .def_ro("reads_passing", &SingleEnzymeQCResult::reads_passing)
        .def_ro("pass_rate", &SingleEnzymeQCResult::pass_rate)
        .def_ro("est_total_fragments", &SingleEnzymeQCResult::est_total_fragments)
        .def_ro("est_total_pairs", &SingleEnzymeQCResult::est_total_pairs)
        .def_ro("avg_fragments_per_read", &SingleEnzymeQCResult::avg_fragments_per_read)
        .def_ro("avg_pairs_per_read", &SingleEnzymeQCResult::avg_pairs_per_read)
        .def_ro("frag_size_mean", &SingleEnzymeQCResult::frag_size_mean)
        .def_ro("frag_size_median", &SingleEnzymeQCResult::frag_size_median)
        .def_ro("frag_size_min", &SingleEnzymeQCResult::frag_size_min)
        .def_ro("frag_size_max", &SingleEnzymeQCResult::frag_size_max)
        .def_ro("read_length_hist_bins", &SingleEnzymeQCResult::read_length_hist_bins)
        .def_ro("read_length_hist_counts", &SingleEnzymeQCResult::read_length_hist_counts)
        .def_ro("sites_hist_bins", &SingleEnzymeQCResult::sites_hist_bins)
        .def_ro("sites_hist_counts", &SingleEnzymeQCResult::sites_hist_counts)
        .def_ro("frag_size_hist_bins", &SingleEnzymeQCResult::frag_size_hist_bins)
        .def_ro("frag_size_hist_counts", &SingleEnzymeQCResult::frag_size_hist_counts);

    // QC function
    m.def("run_qc_analysis_custom", &run_qc_analysis_custom,
          nb::arg("input_path"),
          nb::arg("site"),
          nb::arg("cut_offset"),
          nb::arg("num_reads"),
          nb::arg("min_sites") = 2,
          "Run QC analysis for a custom enzyme site.\n"
          "site: recognition sequence (supports IUPAC degenerate bases)\n"
          "cut_offset: position within site where cut occurs (0 to len(site))\n"
          "Returns SingleEnzymeQCResult with all metrics.");

    // Processing functions
    m.def("process_reads", &process_reads,
          nb::arg("input_path"),
          nb::arg("output_r1"),
          nb::arg("output_r2"),
          nb::arg("enzyme"),
          nb::arg("min_fragments") = 3,
          nb::arg("min_frag_len") = 20,
          nb::arg("strip_overhang") = true,
          nb::arg("gzip_output") = false,
          nb::arg("fast_mode") = false,
          "Process FASTQ or BAM file, generating ALL pairwise contacts (n choose 2).");

    m.def("process_reads_custom", &process_reads_custom,
          nb::arg("input_path"),
          nb::arg("output_r1"),
          nb::arg("output_r2"),
          nb::arg("site"),
          nb::arg("cut_offset"),
          nb::arg("min_fragments") = 3,
          nb::arg("min_frag_len") = 20,
          nb::arg("strip_overhang") = true,
          nb::arg("gzip_output") = false,
          nb::arg("fast_mode") = false,
          "Process FASTQ or BAM file with custom enzyme site.");

    // Enzyme utilities
    m.def("list_enzymes", &cifi::list_enzymes, "Get list of available enzyme names");

    m.def("get_enzyme_info", [](const std::string& name) -> std::pair<std::string, int> {
        auto enzyme_opt = cifi::get_enzyme(name);
        if (!enzyme_opt) {
            throw std::runtime_error("Unknown enzyme: " + name);
        }
        return {enzyme_opt->site, enzyme_opt->cut_offset};
    }, nb::arg("name"), "Get (recognition_site, cut_offset) for enzyme");

    m.def("is_bam_file", &is_bam_file, nb::arg("path"),
          "Check if file is BAM/SAM/CRAM based on extension");

    // IUPAC degenerate base utilities
    m.def("find_all_degenerate", &cifi::find_all_degenerate,
          nb::arg("text"), nb::arg("pattern"),
          "Find all occurrences of pattern in text, supporting IUPAC degenerate bases.");

    m.def("has_degenerate_bases", &cifi::has_degenerate_bases,
          nb::arg("pattern"),
          "Check if pattern contains IUPAC degenerate base codes.");

    m.def("revcomp", &cifi::revcomp, nb::arg("seq"),
          "Reverse complement a DNA sequence.");

    m.def("revcomp_degenerate", &cifi::revcomp_degenerate, nb::arg("seq"),
          "Reverse complement a degenerate DNA sequence (handles IUPAC codes).");

    // BAM filter
    nb::class_<cifi::FilterResult>(m, "FilterResult")
        .def_ro("total_reads", &cifi::FilterResult::total_reads)
        .def_ro("total_pairs", &cifi::FilterResult::total_pairs)
        .def_ro("passed_pairs", &cifi::FilterResult::passed_pairs)
        .def_ro("failed_unpaired", &cifi::FilterResult::failed_unpaired)
        .def_ro("failed_unmapped", &cifi::FilterResult::failed_unmapped)
        .def_ro("failed_mapq", &cifi::FilterResult::failed_mapq)
        .def_ro("failed_mate_not_found", &cifi::FilterResult::failed_mate_not_found);

    m.def("filter_bam", &cifi::filter_bam,
          nb::arg("input_path"),
          nb::arg("output_path"),
          nb::arg("mapq_threshold"),
          nb::arg("threads") = 4,
          "Filter BAM by MAPQ quality for paired reads");
}
