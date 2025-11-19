#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright 2015 Lina Sieverling, Philip Ginsbach, Lars Feuerbach
# C++ wrapper integration

# This file is part of TelomereHunter.

# TelomereHunter is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# TelomereHunter is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with TelomereHunter.  If not, see <http://www.gnu.org/licenses/>.

"""
Python wrapper for C++ telomere filter module
Calls the compiled C++ binary for faster processing
"""

import os
import sys
import subprocess

# Path to C++ binary
CPP_FILTER_BIN = "/igor-shared/TELOMER/telomerehunter_cpp/bin/telomere_filter"

# Fallback to Python implementation if C++ binary not found
try:
    import filter_telomere_reads as python_filter
    HAS_PYTHON_FALLBACK = True
except ImportError:
    HAS_PYTHON_FALLBACK = False


def filter_telomere_reads(bam_file, band_file, out_dir, pid, sample,
                         repeat_threshold_calc, repeat_threshold_set,
                         mapq_threshold, repeats, consecutive_flag,
                         remove_duplicates):
    """
    Wrapper function that calls C++ telomere filter if available,
    otherwise falls back to Python implementation

    Parameters match the original Python implementation
    """

    # Check if C++ binary exists
    if not os.path.exists(CPP_FILTER_BIN):
        print "Warning: C++ filter binary not found at: " + CPP_FILTER_BIN
        if HAS_PYTHON_FALLBACK:
            print "Falling back to Python implementation..."
            return python_filter.filter_telomere_reads(
                bam_file=bam_file,
                band_file=band_file,
                out_dir=out_dir,
                pid=pid,
                sample=sample,
                repeat_threshold_calc=repeat_threshold_calc,
                repeat_threshold_set=repeat_threshold_set,
                mapq_threshold=mapq_threshold,
                repeats=repeats,
                consecutive_flag=consecutive_flag,
                remove_duplicates=remove_duplicates
            )
        else:
            print "Error: No fallback Python implementation available"
            sys.exit(1)

    # Build command for C++ binary
    cmd = [
        CPP_FILTER_BIN,
        "-i", bam_file,
        "-b", band_file,
        "-o", out_dir,
        "-p", pid,
        "-s", sample,
        "-t", str(repeat_threshold_set),
        "-m", str(mapq_threshold),
        "-j", "20"  # Use 20 threads for C++ version
    ]

    # Add repeats
    if repeats:
        cmd.extend(["-r", ",".join(repeats)])

    # Add consecutive flag
    if consecutive_flag:
        cmd.append("-c")

    # Add remove duplicates flag
    if remove_duplicates:
        cmd.append("-d")

    # Print command for debugging
    print "Running C++ filter: " + " ".join(cmd)

    # Execute C++ binary
    try:
        # Run with suppressed progress output (only show errors)
        with open(os.devnull, 'w') as devnull:
            process = subprocess.Popen(
                cmd,
                stdout=devnull,
                stderr=subprocess.PIPE,
                bufsize=1,
                universal_newlines=True
            )

            # Only print error messages, not progress
            for line in iter(process.stderr.readline, ''):
                if line and ('Error' in line or 'Warning' in line or 'Failed' in line):
                    sys.stderr.write(line)
                    sys.stderr.flush()

            # Wait for completion
            return_code = process.wait()

        if return_code != 0:
            print "Error: C++ filter exited with code: " + str(return_code)
            if HAS_PYTHON_FALLBACK:
                print "Falling back to Python implementation..."
                return python_filter.filter_telomere_reads(
                    bam_file=bam_file,
                    band_file=band_file,
                    out_dir=out_dir,
                    pid=pid,
                    sample=sample,
                    repeat_threshold_calc=repeat_threshold_calc,
                    repeat_threshold_set=repeat_threshold_set,
                    mapq_threshold=mapq_threshold,
                    repeats=repeats,
                    consecutive_flag=consecutive_flag,
                    remove_duplicates=remove_duplicates
                )
            else:
                sys.exit(1)

        print "C++ filter completed successfully"

    except Exception as e:
        print "Error running C++ filter: " + str(e)
        if HAS_PYTHON_FALLBACK:
            print "Falling back to Python implementation..."
            return python_filter.filter_telomere_reads(
                bam_file=bam_file,
                band_file=band_file,
                out_dir=out_dir,
                pid=pid,
                sample=sample,
                repeat_threshold_calc=repeat_threshold_calc,
                repeat_threshold_set=repeat_threshold_set,
                mapq_threshold=mapq_threshold,
                repeats=repeats,
                consecutive_flag=consecutive_flag,
                remove_duplicates=remove_duplicates
            )
        else:
            print "Error: No fallback available"
            sys.exit(1)
