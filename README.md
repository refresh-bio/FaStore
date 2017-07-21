# FaStore
FaStore is a high-performance short FASTQ sequencing reads compressor.

The compression is a multi-step process, hence the compressor currently consist of 3 tools:

* `fastore_bin` - performing a DNA reads distribution into bins,
* `fastore_rebin` - performing a DNA reads re-distribution into bins with further clusterization,
* `fastore_pack` - performing compression of the DNA reads stored in bins.

However, for ease of use automated scripts to perform compression and decompression are provided in `scripts` directory. 


# Usage

FaStore offers a variety of different compression configurations. Hence, for an easier selection, we created 4 profiles, namely _lossless_, _reduced_, _lossy_ and _max_. For each of the profiles, a script `fastore_compress_*.sh` is provided.

For example, to compress a pair of FASTQ files `IN_1.fastq` and `IN_2.fastq` with reads represented in pared-end mode and using `8` processing threads type:

    ./compress_*_pe.sh IN_1.fastq IN_2.fastq COMP 8

the compressed files will be stored as `COMP.cmeta' and `COMP.cdata' files.


To decompress the archives generated in any of the above mentioned profile and using `8` processing threads type:

    ./decompress_pe.sh COMP OUT 8

the decompressed files will be stored as `OUT_1.fastq` and `OUT_2.fastq` files.



# Building

## Prerequisites

FaStore currently provides Makefiles for building on Linux and Mac OSX platforms.
However it should also be able to be compiled on Windows platform. 
The only one prerequisite is the availability of the _zlib_ library in the system.


## Compiling

FaStore binaries are by default compiled using g++ version >= 4.8 (reguired for _C++11_ threading support). 
To compile FaStore using _g++_ >= 4.8 with _C++11_ standard and dynamic linking, use the default `Makefile` file and in the main directory type:
    
    make

Alternatively, to compile using clang, invoke `make` using `Makefile.clang` file:

	make -f Makefile.clang

The resulting `fastore_bin`, `fastore_rebin` and `fastore_pack` binaries will be placed in `bin` subdirectory.


However, to compile each subprogram separately, use the makefile files provided in each of subprograms directory.

