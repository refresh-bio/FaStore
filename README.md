# FaStore
FaStore is a set of 3 tools designed for an effective and high-performance compression of short FASTQ sequencing reads, which consist of:
* _fastore\_bin_ - performing a DNA reads distribution into bins,
* _fastore\_rebin_ - performing a DNA reads re-distribution into bins with further clusterization,
* _fastore\_pack_ - performing compression of the DNA reads stored in bins.


# Building

## Prerequisites

FaStore currently provides Makefiles for building on Linux and Mac OSX platforms, however it should also be able to be compiled on Windows platform. The only one prequisite is the availability of the _zlib_ library in the system.


## Compiling

FaStore binaries are by default compiled using g++ version >= 4.8 (reguired for _C++11_ threading support). To compile FaStore using _g++_ >= 4.8 with _C++11_ standard and dynamic linking, in the main directory type:
    
    make

Alternatively, to compile using clang, invoke `make` using `Makefile.clang` file:

	make -f Makefile.clang

The resulting _fastore\_bin_, _fastore\_rebin_ and _fastore\_pack_ binaries will be placed in _bin_ subdirectory.


However, to compile each subprogram separately, use the makefile files provided in each of subprograms directory.


# Usage

FASTQ reads compression using FaStore is a 2-3 stage process, consisting of running _fastore\_bin_ , optionally _fastore\_pack_ and _fastore\_pack_ subprograms chained together. However, to decompress the DNA stream, only running _fastore\_pack_ is needed.

## _fastore\_bin_

...

## _fastore\_rebin_

...

## _fastore\_pack_

...


# Examples

...
