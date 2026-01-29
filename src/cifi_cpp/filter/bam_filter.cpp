#include "bam_filter.hpp"

#include <iostream>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>

namespace cifi {

// Define a struct to hold read information
struct ReadData
{
    bam1_t *bam_record;
    uint32_t flag;
    uint8_t mapq;
    hts_pos_t pos;
    hts_pos_t mpos;
    int32_t tid;
    int32_t mtid;
};

// Hash function for read name prefix
struct PrefixHash
{
    size_t operator()(const std::string &key) const
    {
        return std::hash<std::string>()(key);
    }
};

class BamFilter
{
public:
    BamFilter(const char *input_bam, const char *output_bam, int mapq_cutoff, int n_threads);
    FilterResult process();

private:
    const char *input_bam_;
    const char *output_bam_;
    int mapq_cutoff_;
    int n_threads_;

    // Filtration statistics
    std::atomic<int64_t> total_reads_;
    std::atomic<int64_t> total_pairs_;
    std::atomic<int64_t> passed_pairs_;
    std::atomic<int64_t> failed_unmapped_;
    std::atomic<int64_t> failed_not_paired_;
    std::atomic<int64_t> failed_mapq_;
    std::atomic<int64_t> failed_mate_not_found_;

    // Mutex for hash table access
    std::mutex mtx_;

    // Hash table to store unmatched reads
    std::unordered_map<std::string, std::pair<ReadData, bool>, PrefixHash> read_map_;

    // Maximum number of reads to keep in memory
    size_t max_reads_in_memory_;

    // HTSlib file pointers and header
    htsFile *in_;
    htsFile *out_;
    bam_hdr_t *header_;

    void process_read(bam1_t *read);
    std::string extract_common_prefix(const std::string &qname);
};

BamFilter::BamFilter(const char *input_bam, const char *output_bam, int mapq_cutoff, int n_threads)
    : input_bam_(input_bam), output_bam_(output_bam), mapq_cutoff_(mapq_cutoff), n_threads_(n_threads),
      total_reads_(0), total_pairs_(0), passed_pairs_(0), failed_unmapped_(0), failed_not_paired_(0),
      failed_mapq_(0), failed_mate_not_found_(0),
      max_reads_in_memory_(1000000000), in_(nullptr), out_(nullptr), header_(nullptr)
{
}

FilterResult BamFilter::process()
{
    // Open input BAM file
    in_ = hts_open(input_bam_, "r");
    if (in_ == NULL)
    {
        throw std::runtime_error("Error opening input BAM file");
    }
    // Set threads for reading
    hts_set_threads(in_, n_threads_);

    // Read header from input BAM file
    header_ = sam_hdr_read(in_);
    if (header_ == NULL)
    {
        hts_close(in_);
        throw std::runtime_error("Error reading header from input BAM file");
    }

    // Open output BAM file
    out_ = hts_open(output_bam_, "wb");
    if (out_ == NULL)
    {
        bam_hdr_destroy(header_);
        hts_close(in_);
        throw std::runtime_error("Error opening output BAM file");
    }
    // Set threads for writing
    hts_set_threads(out_, n_threads_);

    // Write header to output BAM file
    if (sam_hdr_write(out_, header_) < 0)
    {
        hts_close(out_);
        bam_hdr_destroy(header_);
        hts_close(in_);
        throw std::runtime_error("Error writing header to output BAM file");
    }

    bam1_t *read = bam_init1();

    while (sam_read1(in_, header_, read) >= 0)
    {
        total_reads_++;
        process_read(read);
        // Reuse the bam1_t struct
        bam_destroy1(read);
        read = bam_init1();
    }

    // Clean up unmatched reads
    for (auto &pair : read_map_)
    {
        bam_destroy1(pair.second.first.bam_record);
        failed_mate_not_found_++;
    }
    read_map_.clear();

    bam_destroy1(read);
    bam_hdr_destroy(header_);
    hts_close(in_);
    hts_close(out_);

    // Build result
    FilterResult result;
    result.total_reads = total_reads_;
    result.total_pairs = total_pairs_;
    result.passed_pairs = passed_pairs_;
    result.failed_unpaired = failed_not_paired_;
    result.failed_unmapped = failed_unmapped_;
    result.failed_mapq = failed_mapq_;
    result.failed_mate_not_found = failed_mate_not_found_;

    return result;
}

std::string BamFilter::extract_common_prefix(const std::string &qname)
{
    // Extract the common prefix up to "ccs:"
    size_t pos = qname.find("ccs:");
    if (pos != std::string::npos)
    {
        return qname.substr(0, pos + 4); // Include "ccs:"
    }
    else
    {
        // Handle cases where "ccs:" is not found
        return qname;
    }
}

void BamFilter::process_read(bam1_t *read)
{
    // Get read name
    std::string qname = bam_get_qname(read);
    std::string prefix = extract_common_prefix(qname);

    // Retrieve flags
    uint32_t flag = read->core.flag;

    // Check if read is paired
    bool is_paired = (flag & BAM_FPAIRED) != 0;

    // Skip reads that are not paired
    if (!is_paired)
    {
        failed_not_paired_++;
        return;
    }

    // Determine if read is first or second in pair
    bool is_first_in_pair = (flag & BAM_FREAD1) != 0;

    // Lock the mutex to modify the hash table
    std::unique_lock<std::mutex> lock(mtx_);

    // Check if the mate is already in the map
    auto it = read_map_.find(prefix);
    if (it != read_map_.end())
    {
        // Mate found, process the pair
        bam1_t *mate = it->second.first.bam_record;
        bool mate_is_first = it->second.second;

        // Ensure that one read is first in pair and the other is second
        if (is_first_in_pair == mate_is_first)
        {
            // Both reads have the same first/second designation, cannot be mates
            failed_mate_not_found_++;
            bam_destroy1(it->second.first.bam_record);
            read_map_.erase(it);
            return;
        }

        total_pairs_++;

        // Retrieve flags
        uint32_t flag1 = read->core.flag;
        uint32_t flag2 = mate->core.flag;

        // Check if reads are unmapped
        bool unmapped1 = (flag1 & BAM_FUNMAP) != 0;
        bool unmapped2 = (flag2 & BAM_FUNMAP) != 0;

        // Get MAPQ scores
        uint8_t mapq1 = read->core.qual;
        uint8_t mapq2 = mate->core.qual;

        // Apply filters
        if (unmapped1 || unmapped2)
        {
            failed_unmapped_++;
        }
        else if (mapq1 < mapq_cutoff_ || mapq2 < mapq_cutoff_)
        {
            failed_mapq_++;
        }
        else
        {
            // Passed all filters, write reads to output BAM file
            if (sam_write1(out_, header_, read) < 0)
            {
                throw std::runtime_error("Error writing to output BAM file");
            }
            if (sam_write1(out_, header_, mate) < 0)
            {
                throw std::runtime_error("Error writing to output BAM file");
            }
            passed_pairs_++;
        }

        // Remove the mate from the map and free memory
        bam_destroy1(mate);
        read_map_.erase(it);
    }
    else
    {
        // Mate not found, store the read
        bam1_t *read_copy = bam_dup1(read);
        read_map_[prefix] = {{read_copy, flag, read->core.qual, read->core.pos, read->core.mpos, read->core.tid, read->core.mtid}, is_first_in_pair};
    }

    // If the map size exceeds the limit, clear it
    if (read_map_.size() > max_reads_in_memory_)
    {
        std::cerr << "Warning: Hash table size exceeded limit. Clearing unmatched reads to free memory.\n";
        // Free memory for unmatched reads
        for (auto &pair : read_map_)
        {
            bam_destroy1(pair.second.first.bam_record);
            failed_mate_not_found_++;
        }
        read_map_.clear();
    }
}

// Wrapper function for Python bindings
FilterResult filter_bam(
    const std::string& input_path,
    const std::string& output_path,
    int mapq_threshold,
    int threads)
{
    BamFilter filter(input_path.c_str(), output_path.c_str(), mapq_threshold, threads);
    return filter.process();
}

}  // namespace cifi
