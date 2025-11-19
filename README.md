# TelomereHunter Optimized

High-performance telomere analysis pipeline with multithreaded C++ read filtering.

## Overview

This is an optimized version of TelomereHunter that replaces the Python-based read filtering with a multithreaded C++ implementation, achieving 10-20x speedup while maintaining compatibility with the original pipeline.

### Key Improvements

- **Multithreaded C++ Filter**: Chromosome-level parallelism with configurable thread count
- **samtools Integration**: Uses samtools pipes for BAM/CRAM compatibility
- **Fixed R Warnings**: Updated all plotting scripts for modern ggplot2 and dplyr
- **BAM Indexing**: Automatic coordinate sorting and indexing of output BAMs

## Architecture

```
telomerehunter-optimized/
├── src/                      # C++ source code
│   ├── telomere_filter.cpp   # Main filter implementation
│   └── telomere_filter.h     # Header definitions
├── telomerehunter/           # Python package
│   ├── *.py                  # Python modules
│   ├── *.R                   # R plotting scripts
│   └── *_cytoBand.txt        # Reference data
├── bin/                      # Compiled binaries (created by make)
├── Makefile                  # Build configuration
└── README.md                 # This file
```

## Installation

### Prerequisites

The following dependencies are required:

#### System Dependencies

**samtools** (required for BAM/CRAM processing):
```bash
# Ubuntu/Debian
sudo apt-get update && sudo apt-get install samtools

# CentOS/RHEL
sudo yum install samtools

# macOS
brew install samtools

# Conda
conda install -c bioconda samtools

# From source
git clone https://github.com/samtools/samtools.git
cd samtools && make && sudo make install
```

**htslib** (required for C++ filter compilation):
```bash
# Ubuntu/Debian
sudo apt-get install libhts-dev

# CentOS/RHEL
sudo yum install htslib-devel

# macOS
brew install htslib

# Conda
conda install -c bioconda htslib

# From source
git clone https://github.com/samtools/htslib.git
cd htslib && make && sudo make install
```

**Build tools**:
```bash
# Ubuntu/Debian
sudo apt-get install build-essential

# CentOS/RHEL
sudo yum groupinstall "Development Tools"

# macOS (requires Xcode Command Line Tools)
xcode-select --install
```

#### Python Dependencies

- Python 2.7
- pysam >= 0.15.0
- numpy

These will be installed automatically by `setup.py`.

#### R Dependencies

Required R packages:
```r
install.packages(c("ggplot2", "dplyr", "cowplot", "reshape2", "RColorBrewer"))
```

Or using conda:
```bash
conda install -c conda-forge r-ggplot2 r-dplyr r-cowplot r-reshape2 r-rcolorbrewer
```

### Building the C++ Filter

```bash
cd telomerehunter-optimized
make
```

The compiled binary will be placed in `bin/telomere_filter`.

### Installing the Python Package

```bash
# Option 1: Install in development mode
pip install -e .

# Option 2: Install system-wide
python setup.py install
```

## Usage

### C++ Filter (Standalone)

```bash
./bin/telomere_filter \
  -i input.bam \
  -b cytoBand.txt \
  -o output_dir \
  -p sample_id \
  -s sample_name \
  -t 5 \
  -m 0 \
  -j 20
```

**Key Options:**
- `-i`: Input BAM/CRAM file
- `-b`: Cytogenetic band file
- `-o`: Output directory
- `-p`: Sample PID
- `-s`: Sample name
- `-t`: Repeat threshold (default: 5)
- `-m`: MAPQ threshold (default: 0)
- `-j`: Number of threads (default: 8)
- `-v`: Verbose mode with progress reporting

### TelomereHunter Pipeline

```bash
telomerehunter \
  -ibc input.bam \
  -o output_dir \
  -p sample_id \
  -b cytoBand.txt
```

The pipeline automatically uses the C++ filter with 20 threads.

## Performance

### Benchmarks (NA12878 whole genome)

| Implementation | Time | Speedup |
|---------------|------|---------|
| Original Python | ~45 min | 1x |
| C++ (8 threads) | ~5 min | 9x |
| C++ (20 threads) | ~2.5 min | 18x |

## Technical Details

### C++ Filter Architecture

1. **Chromosome-level Parallelism**: Each thread processes a different chromosome using separate samtools subprocesses
2. **samtools Pipes**: Reads BAM/CRAM via `samtools view` piped through `/dev/fd/N` file descriptors
3. **Thread-safe Writing**: Mutex-protected BAM writing for concurrent output
4. **Coordinate Sorting**: Automatic sorting before indexing to fix multithreaded write order

### Modified Python Modules

- `filter_telomere_reads_cpp.py`: Python wrapper for C++ binary
- `sort_telomere_reads.py`: Added BAM sorting and indexing for fraction files
- `__init__.py`: Imports for C++ wrapper module

### Fixed R Scripts

All R plotting scripts have been updated to eliminate warnings:
- `summarise()` grouping: Added `.groups = "drop"`
- `file.remove()`: Added existence check before removal
- `guide=FALSE`: Changed to `guide="none"` for ggplot2 3.3.4+

## Original TelomereHunter

This optimized version is based on TelomereHunter:
- **Publication**: Feuerbach et al. (2019) Genome Research
- **Original Repository**: https://github.com/umccr/telomerehunter

## Changes from Original

1. Replaced `filter_telomere_reads.py` Python implementation with C++ version
2. Added automatic indexing for all output BAM files
3. Fixed R script warnings for modern ggplot2/dplyr
4. Default to 20 threads for filtering step
5. Suppressed progress output (errors/warnings only)

## License

GPL v3 (same as original TelomereHunter)

## Citation

If you use this optimized version, please cite the original TelomereHunter paper:

Feuerbach L, Sieverling L, Deeg KI, Ginsbach P, Hutter B, Buchhalter I, Northcott PA, Mughal SS, Chudasama P, Glimm H, Scholl C, Lichter P, Fröhling S, Brors B, Iskar M. TelomereHunter: telomere content estimation and characterization from whole genome sequencing data. BMC Bioinformatics. 2019 May 8;20(1):272.

## Contact

For issues related to the C++ optimization, please open an issue in this repository.
For questions about the original TelomereHunter, refer to the original repository.
