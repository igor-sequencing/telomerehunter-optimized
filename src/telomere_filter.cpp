/*
 * telomere_filter.cpp
 *
 * TelomereHunter filter module implementation
 * Multithreaded BAM/CRAM filtering for telomeric reads
 */

#include "telomere_filter.h"
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <future>
#include <cmath>
#include <regex>
#include <getopt.h>
#include <chrono>
#include <iomanip>
#include <atomic>

// Global mutex for thread-safe output
static std::mutex output_mutex;
static std::mutex bam_write_mutex;

// Global output BAM file
static samFile* global_filtered_bam = nullptr;
static sam_hdr_t* global_bam_header = nullptr;

// Global progress tracking
static std::atomic<uint64_t> global_processed_reads(0);
static std::atomic<uint64_t> global_filtered_reads(0);
static std::atomic<int> global_completed_chromosomes(0);
static int global_total_chromosomes = 0;
static auto global_start_time = std::chrono::steady_clock::now();
static bool global_verbose = false;

// Get reverse complement of DNA sequence
std::string getReverseComplement(const std::string& sequence) {
    std::string temp = sequence;

    // Replace A->1, C->2, G->3, T->4
    for (char& c : temp) {
        switch(c) {
            case 'A': c = '1'; break;
            case 'C': c = '2'; break;
            case 'G': c = '3'; break;
            case 'T': c = '4'; break;
        }
    }

    // Replace 1->T, 2->G, 3->C, 4->A
    for (char& c : temp) {
        switch(c) {
            case '1': c = 'T'; break;
            case '2': c = 'G'; break;
            case '3': c = 'C'; break;
            case '4': c = 'A'; break;
        }
    }

    // Reverse the string
    std::reverse(temp.begin(), temp.end());
    return temp;
}
// Format time duration for display
std::string formatDuration(std::chrono::seconds duration) {
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration % std::chrono::hours(1));
    auto seconds = duration % std::chrono::minutes(1);

    std::ostringstream oss;
    if (hours.count() > 0) {
        oss << hours.count() << "h " << minutes.count() << "m " << seconds.count() << "s";
    } else if (minutes.count() > 0) {
        oss << minutes.count() << "m " << seconds.count() << "s";
    } else {
        oss << seconds.count() << "s";
    }
    return oss.str();
}

// Print overall progress update
void printOverallProgress() {
    if (!global_verbose) return;

    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - global_start_time);

    uint64_t processed = global_processed_reads.load();
    uint64_t filtered = global_filtered_reads.load();
    int completed = global_completed_chromosomes.load();

    double reads_per_sec = processed / (elapsed.count() > 0 ? elapsed.count() : 1);
    double filter_rate = processed > 0 ? (100.0 * filtered / processed) : 0.0;

    std::lock_guard<std::mutex> lock(output_mutex);
    std::cerr << "\n=== Progress Update ===\n";
    std::cerr << "Elapsed time: " << formatDuration(elapsed) << "\n";
    std::cerr << "Chromosomes completed: " << completed << "/" << global_total_chromosomes << "\n";
    std::cerr << "Total reads processed: " << processed << " ("
              << std::fixed << std::setprecision(1) << (reads_per_sec / 1000000.0)
              << "M reads/sec)\n";
    std::cerr << "Telomeric reads found: " << filtered << " ("
              << std::fixed << std::setprecision(2) << filter_rate << "%)\n";
    std::cerr << "======================\n\n";
}



// Calculate GC content percentage (0-100), excluding N bases
int calculateGCContent(const std::string& sequence, float& n_fraction) {
    if (sequence.empty()) {
        n_fraction = 0.0f;
        return 50;
    }

    int gc_count = 0;
    int n_count = 0;

    for (char c : sequence) {
        if (c == 'G' || c == 'C') {
            gc_count++;
        } else if (c == 'N') {
            n_count++;
        }
    }

    int read_length = sequence.length();
    n_fraction = static_cast<float>(n_count) / read_length;

    // Skip if more than 20% N bases
    if (n_fraction > 0.2f) {
        return -1;
    }

    int valid_length = read_length - n_count;
    if (valid_length == 0) {
        return 50;
    }

    return static_cast<int>(round(100.0f * gc_count / valid_length));
}

// Count telomere repeat occurrences
int countTelomereRepeats(const std::string& sequence,
                         const std::vector<std::string>& patterns_forward,
                         const std::vector<std::string>& patterns_reverse,
                         bool consecutive) {
    int max_count = 0;

    if (consecutive) {
        // Search for consecutive repeats using regex
        for (size_t i = 0; i < patterns_forward.size(); i++) {
            const std::string& fwd = patterns_forward[i];
            const std::string& rev = patterns_reverse[i];

            // Create regex for consecutive pattern
            std::string fwd_pattern = "(" + fwd + ")+";
            std::string rev_pattern = "(" + rev + ")+";

            try {
                std::regex fwd_regex(fwd_pattern);
                std::regex rev_regex(rev_pattern);

                std::smatch match;
                std::string::const_iterator searchStart(sequence.cbegin());

                // Count forward pattern
                int fwd_count = 0;
                while (std::regex_search(searchStart, sequence.cend(), match, fwd_regex)) {
                    int match_len = match[0].length();
                    fwd_count = std::max(fwd_count, match_len / static_cast<int>(fwd.length()));
                    searchStart = match.suffix().first;
                }

                // Count reverse pattern
                int rev_count = 0;
                searchStart = sequence.cbegin();
                while (std::regex_search(searchStart, sequence.cend(), match, rev_regex)) {
                    int match_len = match[0].length();
                    rev_count = std::max(rev_count, match_len / static_cast<int>(rev.length()));
                    searchStart = match.suffix().first;
                }

                max_count = std::max(max_count, std::max(fwd_count, rev_count));
            } catch (std::regex_error& e) {
                std::cerr << "Regex error: " << e.what() << std::endl;
            }
        }
    } else {
        // Non-consecutive: count all occurrences
        for (size_t i = 0; i < patterns_forward.size(); i++) {
            const std::string& fwd = patterns_forward[i];
            const std::string& rev = patterns_reverse[i];

            int fwd_count = 0;
            size_t pos = 0;
            while ((pos = sequence.find(fwd, pos)) != std::string::npos) {
                fwd_count++;
                pos += fwd.length();
            }

            int rev_count = 0;
            pos = 0;
            while ((pos = sequence.find(rev, pos)) != std::string::npos) {
                rev_count++;
                pos += rev.length();
            }

            max_count = std::max(max_count, std::max(fwd_count, rev_count));
        }
    }

    return max_count;
}

// Find which chromosome band a position belongs to
std::string findChromosomeBand(const std::vector<ChromosomeBand>& bands,
                               int64_t position) {
    for (const auto& band : bands) {
        if (position <= static_cast<int64_t>(band.end)) {
            return band.band_name;
        }
    }
    return bands.empty() ? "unknown" : bands.back().band_name;
}

// Load chromosome banding file
std::map<std::string, std::vector<ChromosomeBand>>
loadBandingFile(const std::string& band_file, const std::string& chr_prefix) {
    std::map<std::string, std::vector<ChromosomeBand>> bands;

    // Initialize chromosomes
    std::vector<std::string> chromosomes;
    for (int i = 1; i <= 22; i++) {
        chromosomes.push_back(std::to_string(i));
    }
    chromosomes.push_back("X");
    chromosomes.push_back("Y");
    chromosomes.push_back("unmapped");

    for (const auto& chr : chromosomes) {
        bands[chr] = std::vector<ChromosomeBand>();
    }

    // Read banding file
    std::ifstream file(band_file);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open banding file: " << band_file << std::endl;
        return bands;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string chr, start_str, end_str, band_name;

        if (!(iss >> chr >> start_str >> end_str >> band_name)) {
            continue;  // Skip invalid lines
        }

        // Remove "chr" prefix if present
        if (chr.substr(0, 3) == "chr") {
            chr = chr.substr(3);
        }

        uint64_t end = std::stoull(end_str);

        if (bands.find(chr) != bands.end()) {
            bands[chr].emplace_back(band_name, end);
        }
    }

    file.close();

    // Add unmapped band
    bands["unmapped"].emplace_back("unmapped", 0);

    return bands;
}

// Process a single chromosome using samtools view via pipe (thread-safe, works with BAM and CRAM)
FilterResults processChromosome(const std::string& bampath,
                                const std::string& chr_name,
                                int tid,
                                const FilterParams& params,
                                const std::map<std::string, std::vector<ChromosomeBand>>& bands,
                                const std::string& chr_prefix) {
    FilterResults results;

    // Copy band structure for this chromosome
    results.bands = bands;

    {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "[" << chr_name << "] Starting processing via samtools..." << std::endl;
    }

    // Prepare reverse complement patterns
    std::vector<std::string> patterns_reverse;
    for (const auto& pattern : params.repeats) {
        patterns_reverse.push_back(getReverseComplement(pattern));
    }

    // Build samtools command to extract reads for this chromosome
    std::string samtools_cmd;
    if (tid == -1 || chr_name == "*") {
        // Unmapped reads
        samtools_cmd = "samtools view -h \"" + bampath + "\" '*'";
    } else {
        // Specific chromosome
        samtools_cmd = "samtools view -h \"" + bampath + "\" \"" + chr_name + "\"";
    }

    // Open pipe to samtools
    FILE* pipe_fp = popen(samtools_cmd.c_str(), "r");
    if (pipe_fp == NULL) {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "Error: Failed to open samtools pipe for " << chr_name << std::endl;
        return results;
    }

    // Open SAM stream from pipe using /dev/fd/N
    char fd_path[256];
    snprintf(fd_path, sizeof(fd_path), "/dev/fd/%d", fileno(pipe_fp));

    samFile* sam_fp = sam_open(fd_path, "r");
    if (sam_fp == NULL) {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "Error: Failed to open SAM stream for " << chr_name << std::endl;
        pclose(pipe_fp);
        return results;
    }

    // Read header
    sam_hdr_t* header = sam_hdr_read(sam_fp);
    if (header == NULL) {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "Error: Failed to read header for " << chr_name << std::endl;
        sam_close(sam_fp);
        pclose(pipe_fp);
        return results;
    }

    // Allocate BAM alignment structure
    bam1_t* b = bam_init1();
    int nprocessed = 0;
    int nfiltered = 0;

    // Get chromosome name without prefix for band lookup
    std::string chr_for_bands = chr_name;
    if (chr_name == "*") {
        chr_for_bands = "unmapped";
    } else if (chr_for_bands.substr(0, chr_prefix.length()) == chr_prefix) {
        chr_for_bands = chr_for_bands.substr(chr_prefix.length());
    }

    // Process reads from samtools pipe
    int ret;
    while ((ret = sam_read1(sam_fp, header, b)) >= 0) {
        // Skip secondary alignments
        if (b->core.flag & BAM_FSECONDARY) continue;

        // Skip duplicates if requested
        if (params.remove_duplicates && (b->core.flag & BAM_FDUP)) continue;

        // Skip supplementary alignments
        if (b->core.flag & BAM_FSUPPLEMENTARY) continue;

        // Extract sequence
        std::string sequence;
        uint8_t* seq = bam_get_seq(b);
        for (int i = 0; i < b->core.l_qseq; i++) {
            sequence += seq_nt16_str[bam_seqi(seq, i)];
        }

        if (sequence.empty()) continue;

        results.total_reads++;
        nprocessed++;

        // Calculate GC content
        float n_fraction;
        int gc_content = calculateGCContent(sequence, n_fraction);

        if (gc_content >= 0 && n_fraction <= 0.2f) {
            results.gc_content[gc_content]++;
        }

        // Determine chromosome and band for read counting
        std::string read_chr = "unmapped";
        std::string band_name = "unmapped";

        bool is_unmapped = (b->core.flag & BAM_FUNMAP) || (b->core.qual < params.mapq_threshold);

        if (!is_unmapped && tid != -1) {
            read_chr = chr_for_bands;
            int64_t pos = b->core.pos;

            if (results.bands.find(read_chr) != results.bands.end() &&
                !results.bands[read_chr].empty()) {
                band_name = findChromosomeBand(results.bands[read_chr], pos);

                // Increment read count for this band
                for (auto& band : results.bands[read_chr]) {
                    if (band.band_name == band_name) {
                        band.read_count++;
                        break;
                    }
                }
            }
        } else {
            // Unmapped read
            if (results.bands.find("unmapped") != results.bands.end() &&
                !results.bands["unmapped"].empty()) {
                results.bands["unmapped"][0].read_count++;
            }
        }

        // Calculate repeat threshold
        int repeat_threshold = params.repeat_threshold_set;
        if (params.repeat_threshold_calc == -1) {  // heterogeneous
            repeat_threshold = static_cast<int>(round(sequence.length() * params.repeat_threshold_set / 100.0));
        }

        // Count telomere repeats
        int repeat_count = countTelomereRepeats(sequence, params.repeats,
                                                patterns_reverse, params.consecutive_flag);

        // Check if read passes threshold
        if (repeat_count >= repeat_threshold) {
            results.telomeric_reads++;
            results.filtered_reads++;
            nfiltered++;

            // Write to global output BAM (thread-safe)
            {
                std::lock_guard<std::mutex> lock(bam_write_mutex);
                if (global_filtered_bam != nullptr && global_bam_header != nullptr) {
                    if (sam_write1(global_filtered_bam, global_bam_header, b) < 0) {
                        std::cerr << "Warning: Failed to write filtered read" << std::endl;
                    }
                }
            }
        }

        if (params.verbose) {
            if (nprocessed % 500000 == 0 && nprocessed > 0) {
                std::lock_guard<std::mutex> lock(output_mutex);
                std::cerr << "[" << chr_name << "] Progress: " << nprocessed
                          << " reads processed, " << nfiltered << " telomeric ("
                          << std::fixed << std::setprecision(2)
                          << (100.0 * nfiltered / nprocessed) << "%)\n";
            }
        } else {
            if (nprocessed % 1000000 == 0) {
                std::lock_guard<std::mutex> lock(output_mutex);
                std::cerr << "[" << chr_name << "] Processed " << nprocessed
                          << " reads, filtered " << nfiltered << std::endl;
            }
        }
    }

    // Cleanup
    bam_destroy1(b);
    sam_hdr_destroy(header);
    sam_close(sam_fp);  // Close htslib handle first
    pclose(pipe_fp);    // Then close the pipe

    // Update global counters
    global_processed_reads += nprocessed;
    global_filtered_reads += nfiltered;
    global_completed_chromosomes++;

    {
        std::lock_guard<std::mutex> lock(output_mutex);
        std::cerr << "[" << chr_name << "] Completed. Processed " << nprocessed
                  << " reads, filtered " << nfiltered << std::endl;
    }

    if (params.verbose) {
        printOverallProgress();
    }

    return results;
}

// Merge results from multiple chromosomes
void mergeResults(FilterResults& dest, const FilterResults& src) {
    dest.total_reads += src.total_reads;
    dest.telomeric_reads += src.telomeric_reads;
    dest.filtered_reads += src.filtered_reads;

    // Merge GC content
    for (const auto& pair : src.gc_content) {
        dest.gc_content[pair.first] += pair.second;
    }

    // Merge band counts
    for (const auto& chr_pair : src.bands) {
        const std::string& chr = chr_pair.first;
        if (dest.bands.find(chr) == dest.bands.end()) {
            dest.bands[chr] = chr_pair.second;
        } else {
            for (size_t i = 0; i < chr_pair.second.size() && i < dest.bands[chr].size(); i++) {
                dest.bands[chr][i].read_count += chr_pair.second[i].read_count;
            }
        }
    }
}

// Main filtering function
int filterTelomereReads(const FilterParams& params) {
    std::cerr << "Starting telomere read filtering..." << std::endl;
    std::cerr << "Input BAM: " << params.bam_file << std::endl;
    std::cerr << "Output directory: " << params.out_dir << std::endl;
    std::cerr << "Using " << params.num_threads << " threads" << std::endl;

    // Open input BAM to read header
    samFile* sam_fp = sam_open(params.bam_file.c_str(), "r");
    if (sam_fp == NULL) {
        std::cerr << "Error: Could not open BAM file: " << params.bam_file << std::endl;
        return 1;
    }

    sam_hdr_t* header = sam_hdr_read(sam_fp);
    if (header == NULL) {
        std::cerr << "Error: Could not read BAM header" << std::endl;
        sam_close(sam_fp);
        return 1;
    }

    // Detect chromosome prefix
    std::string chr_prefix = "";
    if (header->n_targets > 0) {
        std::string first_chr = header->target_name[0];
        if (first_chr.substr(0, 3) == "chr") {
            chr_prefix = "chr";
        }
    }

    std::cerr << "Detected chromosome prefix: '" << chr_prefix << "'" << std::endl;

    // Load banding file
    auto bands = loadBandingFile(params.band_file, chr_prefix);
    std::cerr << "Loaded " << bands.size() << " chromosome banding regions" << std::endl;

    // Open output BAM file
    std::string filtered_bam_path = params.out_dir + "/" + params.pid + "_filtered.bam";
    global_filtered_bam = sam_open(filtered_bam_path.c_str(), "wb");
    if (global_filtered_bam == NULL) {
        std::cerr << "Error: Could not create output BAM file" << std::endl;
        sam_hdr_destroy(header);
        sam_close(sam_fp);
        return 1;
    }

    global_bam_header = sam_hdr_dup(header);
    if (sam_hdr_write(global_filtered_bam, global_bam_header) < 0) {
        std::cerr << "Error: Could not write BAM header" << std::endl;
        sam_close(global_filtered_bam);
        sam_hdr_destroy(global_bam_header);
        sam_hdr_destroy(header);
        sam_close(sam_fp);
        return 1;
    }

    // Check for index
    hts_idx_t* idx = sam_index_load(sam_fp, params.bam_file.c_str());
    if (idx == NULL) {
        std::cerr << "Error: Index file not found. Please index your BAM file." << std::endl;
        sam_close(global_filtered_bam);
        sam_hdr_destroy(global_bam_header);
        sam_hdr_destroy(header);
        sam_close(sam_fp);
        return 1;
    }

    // Collect chromosomes to process
    std::vector<std::pair<int, std::string>> chromosomes;

    // Add unmapped reads
    chromosomes.push_back({-1, "*"});

    // Add all chromosomes from header
    for (int i = 0; i < header->n_targets; i++) {
        chromosomes.push_back({i, header->target_name[i]});
    }

    hts_idx_destroy(idx);
    sam_hdr_destroy(header);
    sam_close(sam_fp);

    std::cerr << "Processing " << chromosomes.size() << " regions in parallel..." << std::endl;

    // Initialize progress tracking
    global_start_time = std::chrono::steady_clock::now();
    global_processed_reads = 0;
    global_filtered_reads = 0;
    global_completed_chromosomes = 0;
    global_total_chromosomes = chromosomes.size();
    global_verbose = params.verbose;

    if (params.verbose) {
        std::cerr << "\n=== Starting Telomere Filtering ===\n";
        std::cerr << "Total regions to process: " << global_total_chromosomes << "\n";
        std::cerr << "Threads: " << params.num_threads << "\n";
        std::cerr << "Repeat threshold: " << params.repeat_threshold_set << "\n";
        std::cerr << "MAPQ threshold: " << params.mapq_threshold << "\n";
        std::cerr << "==============================\n\n";
    }

    // Process chromosomes in parallel
    FilterResults final_results;
    std::vector<std::future<FilterResults>> futures;

    for (const auto& chr_pair : chromosomes) {
        int tid = chr_pair.first;
        std::string chr_name = chr_pair.second;

        // Limit concurrent threads
        while (futures.size() >= static_cast<size_t>(params.num_threads)) {
            for (auto it = futures.begin(); it != futures.end(); ) {
                if (it->wait_for(std::chrono::milliseconds(10)) == std::future_status::ready) {
                    auto result = it->get();
                    mergeResults(final_results, result);
                    it = futures.erase(it);
                } else {
                    ++it;
                }
            }
            if (futures.size() >= static_cast<size_t>(params.num_threads)) {
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
        }

        // Launch new thread
        futures.push_back(std::async(std::launch::async,
                                     processChromosome,
                                     params.bam_file,
                                     chr_name,
                                     tid,
                                     params,
                                     bands,
                                     chr_prefix));
    }

    // Collect remaining results
    for (auto& future : futures) {
        auto result = future.get();
        mergeResults(final_results, result);
    }

    // Close output BAM
    sam_close(global_filtered_bam);
    sam_hdr_destroy(global_bam_header);
    global_filtered_bam = nullptr;
    global_bam_header = nullptr;

    std::cerr << "Filtered " << final_results.filtered_reads << " telomeric reads out of "
              << final_results.total_reads << " total reads" << std::endl;

    // Write readcount file
    std::string readcount_path = params.out_dir + "/" + params.pid + "_readcount.tsv";
    std::ofstream readcount_file(readcount_path);
    if (readcount_file.is_open()) {
        readcount_file << "chr\tband\treads\n";

        std::vector<std::string> chromosomes_ordered;
        for (int i = 1; i <= 22; i++) {
            chromosomes_ordered.push_back(std::to_string(i));
        }
        chromosomes_ordered.push_back("X");
        chromosomes_ordered.push_back("Y");
        chromosomes_ordered.push_back("unmapped");

        for (const auto& chr : chromosomes_ordered) {
            if (final_results.bands.find(chr) != final_results.bands.end()) {
                for (const auto& band : final_results.bands[chr]) {
                    readcount_file << chr << "\t" << band.band_name << "\t"
                                   << band.read_count << "\n";
                }
            }
        }

        readcount_file.close();
        std::cerr << "Wrote readcount file: " << readcount_path << std::endl;
    }

    // Write GC content file
    std::string gc_content_path = params.out_dir + "/" + params.pid + "_" + params.sample + "_gc_content.tsv";
    std::ofstream gc_content_file(gc_content_path);
    if (gc_content_file.is_open()) {
        gc_content_file << "gc_content_percent\tread_count\n";

        for (int gc = 0; gc <= 100; gc++) {
            gc_content_file << gc << "\t" << final_results.gc_content[gc] << "\n";
        }

        gc_content_file.close();
        std::cerr << "Wrote GC content file: " << gc_content_path << std::endl;
    }

    // Sort filtered BAM by coordinate (required before indexing)
    // The BAM is unsorted because multiple threads wrote to it
    std::cerr << "Sorting filtered BAM by coordinate..." << std::endl;
    std::string sorted_bam_path = params.out_dir + "/" + params.pid + "_filtered_sorted.bam";
    std::string sort_coord_cmd = "samtools sort -@ " + std::to_string(params.num_threads) +
                                  " \"" + filtered_bam_path + "\" -o \"" + sorted_bam_path + "\"";
    int ret = system(sort_coord_cmd.c_str());
    if (ret != 0) {
        std::cerr << "Error: Failed to sort BAM by coordinate" << std::endl;
        return 1;
    }

    // Replace unsorted with sorted
    std::string mv_cmd = "mv \"" + sorted_bam_path + "\" \"" + filtered_bam_path + "\"";
    ret = system(mv_cmd.c_str());
    if (ret != 0) {
        std::cerr << "Error: Failed to replace unsorted BAM" << std::endl;
        return 1;
    }
    std::cerr << "Sorted BAM saved to: " << filtered_bam_path << std::endl;

    // Index the sorted BAM
    std::cerr << "Indexing filtered BAM file..." << std::endl;
    std::string index_cmd = "samtools index \"" + filtered_bam_path + "\"";
    ret = system(index_cmd.c_str());
    if (ret != 0) {
        std::cerr << "Warning: Failed to index filtered BAM" << std::endl;
    } else {
        std::cerr << "Indexed BAM: " << filtered_bam_path << ".bai" << std::endl;
    }

    // Sort by name for downstream processing
    std::cerr << "Sorting filtered BAM by name..." << std::endl;
    std::string name_sorted_path = params.out_dir + "/" + params.pid + "_filtered_name_sorted.bam";
    std::string sort_name_cmd = "samtools sort -n -@ " + std::to_string(params.num_threads) +
                                 " \"" + filtered_bam_path + "\" -o \"" + name_sorted_path + "\"";
    ret = system(sort_name_cmd.c_str());
    if (ret != 0) {
        std::cerr << "Warning: Failed to sort BAM by name" << std::endl;
    } else {
        std::cerr << "Created name-sorted BAM: " << name_sorted_path << std::endl;
    }

    std::cerr << "Filtering complete!" << std::endl;
    return 0;
}

// Main program
static const char* USAGE_MESSAGE =
"Program: telomere_filter\n"
"Usage: telomere_filter [OPTIONS]\n\n"
"Required:\n"
"  -i, --input FILE        Input BAM/CRAM file (indexed)\n"
"  -b, --bands FILE        Chromosome banding file\n"
"  -o, --outdir DIR        Output directory\n"
"  -p, --pid STRING        Sample/patient ID\n"
"  -s, --sample STRING     Sample name (tumor/control)\n\n"
"Optional:\n"
"  -r, --repeats STR       Telomere repeat patterns (comma-separated, default: TTAGGG,TGAGGG,TCAGGG,TTGGGG)\n"
"  -t, --threshold INT     Repeat threshold (default: 6 per 100bp)\n"
"  -m, --mapq INT          Mapping quality threshold (default: 8)\n"
"  -c, --consecutive       Search for consecutive repeats only\n"
"  -d, --remove-duplicates Remove duplicate reads\n"
"  -j, --threads INT       Number of threads (default: 1)\n"
"  -v, --verbose           Enable detailed progress reporting\n"
"  -h, --help              Display this help message\n\n";

int main(int argc, char** argv) {
    FilterParams params;

    static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"bands", required_argument, 0, 'b'},
        {"outdir", required_argument, 0, 'o'},
        {"pid", required_argument, 0, 'p'},
        {"sample", required_argument, 0, 's'},
        {"repeats", required_argument, 0, 'r'},
        {"threshold", required_argument, 0, 't'},
        {"mapq", required_argument, 0, 'm'},
        {"consecutive", no_argument, 0, 'c'},
        {"remove-duplicates", no_argument, 0, 'd'},
        {"threads", required_argument, 0, 'j'},
        {"verbose", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "i:b:o:p:s:r:t:m:cdj:vh", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'i':
                params.bam_file = optarg;
                break;
            case 'b':
                params.band_file = optarg;
                break;
            case 'o':
                params.out_dir = optarg;
                break;
            case 'p':
                params.pid = optarg;
                break;
            case 's':
                params.sample = optarg;
                break;
            case 'r': {
                params.repeats.clear();
                std::string repeats_str = optarg;
                std::istringstream iss(repeats_str);
                std::string repeat;
                while (std::getline(iss, repeat, ',')) {
                    params.repeats.push_back(repeat);
                }
                break;
            }
            case 't':
                params.repeat_threshold_set = std::atoi(optarg);
                break;
            case 'm':
                params.mapq_threshold = std::atoi(optarg);
                break;
            case 'c':
                params.consecutive_flag = true;
                break;
            case 'd':
                params.remove_duplicates = true;
                break;
            case 'j':
                params.num_threads = std::atoi(optarg);
                if (params.num_threads < 1) params.num_threads = 1;
                break;
            case 'v':
                params.verbose = true;
                break;
            case 'h':
                std::cout << USAGE_MESSAGE;
                return 0;
            default:
                std::cerr << USAGE_MESSAGE;
                return 1;
        }
    }

    // Validate required parameters
    if (params.bam_file.empty() || params.band_file.empty() ||
        params.out_dir.empty() || params.pid.empty() || params.sample.empty()) {
        std::cerr << "Error: Missing required parameters\n\n";
        std::cerr << USAGE_MESSAGE;
        return 1;
    }

    // Run filtering
    return filterTelomereReads(params);
}
