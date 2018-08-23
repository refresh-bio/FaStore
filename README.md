# FaStore

## Overview

FaStore is a high-performance short FASTQ sequencing reads compressor.

The compression process happens over multiple steps and the compressor
currently consists of 3 tools:

* `fastore_bin` - distribute DNA reads into bins
* `fastore_rebin` - re-distribute the DNA reads into bins and clusterize further
* `fastore_pack` - compress the DNA reads stored in bins

However, for ease of use, automated scripts are provided in the `scripts`
directory to perform compression and decompression.


## Usage

FaStore offers a variety of different compression configurations. In order to
simplify selection, we created 4 profiles, namely _lossless_, _reduced_,
_lossy_ and _max_. To perform automatic compression and decompression, a pair
of scripts `fastore_compress.sh` and `fastore_decompress.sh` is provided.

For example, to compress a pair of FASTQ files `IN_1.fastq` and `IN_2.fastq` in
the lossless mode with reads represented in pared-end mode and using `8`
processing threads type:

```bash
./fastore_compress.sh --lossless --in IN_1.fastq --pair IN_2.fastq --out COMP --threads 8
```

the compressed files will be stored as `COMP.cmeta` and `COMP.cdata` files.


To decompress the archives generated with any of the above mentioned profiles
and using `8` processing threads type:

```bash
./fastore_decompress.sh --in COMP --out OUT_1.fastq --pair OUT_2.fastq --threads 8
```

the decompressed files will be stored as `OUT_1.fastq` and `OUT_2.fastq` files.


## Building

### Prerequisites

FaStore currently provides Makefiles for building on Linux and Mac OSX
platforms. However, it should also be able to be compiled on Windows platform.

**The only prerequisite is the _zlib_ library.**


### Compiling

FaStore binaries are by default compiled using g++ version >= 4.8 (reguired for
_C++11_ threading support). To compile FaStore using _g++_ >= 4.8 with _C++11_
standard and dynamic linking, use the default `Makefile` file and in the main
directory type:

```bash
make
```

Alternatively, to compile using clang, invoke `make` using `Makefile.clang`
file:

```bash
make -f Makefile.clang
```

The resulting `fastore_bin`, `fastore_rebin` and `fastore_pack` binaries will
be placed in `bin` subdirectory.

However, to compile each subprogram separately, use the makefile files provided
in each of subprograms directory.
