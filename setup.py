#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
TelomereHunter Optimized Setup
High-performance telomere analysis with C++ acceleration
"""

from setuptools import setup, find_packages
from setuptools.command.install import install
import subprocess
import os

class CustomInstall(install):
    """Custom installation to build C++ component"""

    def check_samtools(self):
        """Check if samtools is installed"""
        try:
            subprocess.check_output(['samtools', '--version'], stderr=subprocess.STDOUT)
            print("✓ samtools found")
            return True
        except (OSError, subprocess.CalledProcessError):
            print("\n" + "="*70)
            print("WARNING: samtools is not installed or not in PATH")
            print("="*70)
            print("\nsamtools is required for TelomereHunter to process BAM/CRAM files.")
            print("\nInstallation instructions:")
            print("  Ubuntu/Debian:  sudo apt-get install samtools")
            print("  CentOS/RHEL:    sudo yum install samtools")
            print("  macOS:          brew install samtools")
            print("  Conda:          conda install -c bioconda samtools")
            print("  From source:    https://github.com/samtools/samtools")
            print("\nAfter installing samtools, you can continue using TelomereHunter.")
            print("="*70 + "\n")
            return False

    def check_htslib(self):
        """Check if htslib development files are available"""
        try:
            # Try to compile a simple test
            result = subprocess.check_output(
                ['pkg-config', '--cflags', '--libs', 'htslib'],
                stderr=subprocess.STDOUT
            )
            print("✓ htslib development files found")
            return True
        except (OSError, subprocess.CalledProcessError):
            print("\n" + "="*70)
            print("WARNING: htslib development files not found")
            print("="*70)
            print("\nhtslib is required to build the C++ filter.")
            print("\nInstallation instructions:")
            print("  Ubuntu/Debian:  sudo apt-get install libhts-dev")
            print("  CentOS/RHEL:    sudo yum install htslib-devel")
            print("  macOS:          brew install htslib")
            print("  Conda:          conda install -c bioconda htslib")
            print("  From source:    https://github.com/samtools/htslib")
            print("\nAfter installing htslib, run: python setup.py install")
            print("="*70 + "\n")
            return False

    def run(self):
        # Check dependencies
        print("\nChecking dependencies...")
        samtools_ok = self.check_samtools()
        htslib_ok = self.check_htslib()

        # Build C++ filter
        print("\nBuilding C++ telomere filter...")
        try:
            subprocess.check_call(['make', 'clean'])
            subprocess.check_call(['make'])
            print("✓ C++ filter built successfully\n")
        except subprocess.CalledProcessError as e:
            print("\nWarning: Failed to build C++ filter.")
            print("Error: {}".format(e))
            print("You may need to install dependencies and build manually.")
            print("Run: make clean && make\n")

        # Run standard installation
        install.run(self)

setup(
    name='telomerehunter-optimized',
    version='1.1.0',
    description='High-performance telomere analysis with multithreaded C++ filtering',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Original: Lars Feuerbach et al.; Optimized: C++ acceleration',
    author_email='',
    url='https://github.com/umccr/telomerehunter',
    license='GPLv3',

    packages=['telomerehunter'],
    package_dir={'telomerehunter': 'telomerehunter'},
    package_data={
        'telomerehunter': [
            '*.R',
            '*_cytoBand.txt',
            '*.tsv',
        ],
    },

    scripts=['bin/telomerehunter'],

    data_files=[
        ('bin', ['bin/telomere_filter']),
    ],

    install_requires=[
        'pysam>=0.15.0',
        'numpy',
    ],

    cmdclass={
        'install': CustomInstall,
    },

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: C++',
    ],

    python_requires='>=2.7, <3',

    keywords='bioinformatics telomeres genomics cancer whole-genome-sequencing',

    project_urls={
        'Original TelomereHunter': 'https://github.com/umccr/telomerehunter',
        'Publication': 'https://doi.org/10.1186/s12859-019-2851-0',
    },
)
