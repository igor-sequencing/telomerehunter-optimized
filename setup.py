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
    def run(self):
        # Build C++ filter
        print("Building C++ telomere filter...")
        try:
            subprocess.check_call(['make', 'clean'])
            subprocess.check_call(['make'])
            print("C++ filter built successfully")
        except subprocess.CalledProcessError as e:
            print("Warning: Failed to build C++ filter. You may need to build it manually.")
            print("Error: {}".format(e))

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
