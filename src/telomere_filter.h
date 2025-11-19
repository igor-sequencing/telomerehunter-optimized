/*
 * telomere_filter.h
 *
 * TelomereHunter filter module - C++ implementation
 * Filters telomeric reads from BAM/CRAM files with multithreading support
 *
 * Based on original TelomereHunter Python implementation
 * Copyright 2015 Lina Sieverling, Philip Ginsbach, Lars Feuerbach
 * C++ port with htslib and multithreading
 */

#ifndef TELOMERE_FILTER_H_
#define TELOMERE_FILTER_H_

#include <string>
#include <vector>
#include <map>
#include <mutex>
#include <cstdint>

// Structures

struct ChromosomeBand {
    std::string band_name;
    uint64_t end;
    uint64_t read_count;

    ChromosomeBand() : band_name(""), end(0), read_count(0) {}
    ChromosomeBand(const std::string& name, uint64_t e)
        : band_name(name), end(e), read_count(0) {}
};

struct FilterResults {
    // Per-chromosome band data
    std::map<std::string, std::vector<ChromosomeBand>> bands;

    // GC content distribution (0-100%)
    std::map<int, uint64_t> gc_content;

    // Statistics
    uint64_t total_reads;
    uint64_t telomeric_reads;
    uint64_t filtered_reads;

    FilterResults() : total_reads(0), telomeric_reads(0), filtered_reads(0) {
        // Initialize GC content bins
        for (int i = 0; i <= 100; i++) {
            gc_content[i] = 0;
        }
    }
};

struct FilterParams {
    std::string bam_file;
    std::string band_file;
    std::string out_dir;
    std::string pid;
    std::string sample;

    // Repeat patterns to search for
    std::vector<std::string> repeats;

    // Thresholds
    int repeat_threshold_calc;      // -1 for 'n' (heterogeneous)
    int repeat_threshold_set;       // Set threshold or percentage
    int mapq_threshold;

    // Flags
    bool consecutive_flag;
    bool remove_duplicates;
    bool per_read_length;
    bool verbose;

    // Threading
    int num_threads;

    FilterParams()
        : repeat_threshold_calc(-1),
          repeat_threshold_set(6),
          mapq_threshold(8),
          consecutive_flag(false),
          remove_duplicates(false),
          per_read_length(true),
          verbose(false),
          num_threads(1) {
        // Default repeats
        repeats = {"TTAGGG", "TGAGGG", "TCAGGG", "TTGGGG"};
    }
};

// Function declarations

/**
 * Get reverse complement of DNA sequence
 */
std::string getReverseComplement(const std::string& sequence);

/**
 * Calculate GC content percentage (0-100)
 * Excludes N bases from calculation
 */
int calculateGCContent(const std::string& sequence, float& n_fraction);

/**
 * Count occurrences of telomere patterns in sequence
 * Returns max count of forward or reverse complement
 */
int countTelomereRepeats(const std::string& sequence,
                         const std::vector<std::string>& patterns_forward,
                         const std::vector<std::string>& patterns_reverse,
                         bool consecutive);

/**
 * Find which chromosome band a position belongs to
 */
std::string findChromosomeBand(const std::vector<ChromosomeBand>& bands,
                               int64_t position);

/**
 * Load chromosome banding file
 */
std::map<std::string, std::vector<ChromosomeBand>>
loadBandingFile(const std::string& band_file,
                const std::string& chr_prefix);

/**
 * Process a single chromosome in parallel
 */
FilterResults processChromosome(const std::string& bampath,
                                const std::string& chr_name,
                                int tid,
                                const FilterParams& params,
                                const std::map<std::string, std::vector<ChromosomeBand>>& bands,
                                const std::string& chr_prefix);

/**
 * Main filtering function - processes entire BAM file
 */
int filterTelomereReads(const FilterParams& params);

/**
 * Write output files: filtered BAM, readcount.tsv, gc_content.tsv
 */
bool writeOutputFiles(const FilterParams& params,
                      const FilterResults& results,
                      const std::vector<std::string>& filtered_read_names);

/**
 * Sort and index BAM file
 */
bool sortAndIndexBam(const std::string& bam_file, const std::string& output_prefix);

#endif /* TELOMERE_FILTER_H_ */
