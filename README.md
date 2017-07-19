# FaStore
FaStore is a high-performance short FASTQ sequencing reads compressor.

The compression is a multi-step process, hence the compressor consist of 3 tools:

* `fastore_bin` - performing a DNA reads distribution into bins,
* `fastore_rebin` - performing a DNA reads re-distribution into bins with further clusterization,
* `fastore_pack` - performing compression of the DNA reads stored in bins.


# Building

## Prerequisites

FaStore currently provides Makefiles for building on Linux and Mac OSX platforms.
However it should also be able to be compiled on Windows platform. 
The only one prequisite is the availability of the _zlib_ library in the system.


## Compiling

FaStore binaries are by default compiled using g++ version >= 4.8 (reguired for _C++11_ threading support). 
To compile FaStore using _g++_ >= 4.8 with _C++11_ standard and dynamic linking, use the default `Makefile` file and in the main directory type:
    
    make

Alternatively, to compile using clang, invoke `make` using `Makefile.clang` file:

	make -f Makefile.clang

The resulting `fastore_bin_`, `_fastore_rebin` and `fastore_pack` binaries will be placed in _bin_ subdirectory.


However, to compile each subprogram separately, use the makefile files provided in each of subprograms directory.


# Usage

FASTQ reads compression is a multi-step process, consisting of running `fastore_bin` , optionally `fastore_rebin` and `fastore_pack` subprograms chained together. However, to decompress the DNA stream, only running _fastore\_pack_ is needed.

## _fastore\_bin_

...

## _fastore\_rebin_

...

## _fastore\_pack_

...


# Examples

...
