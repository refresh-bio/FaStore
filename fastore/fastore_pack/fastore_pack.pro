TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

QMAKE_CXXFLAGS += -std=c++11
QMAKE_CXXFLAGS += -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE

LIBS += -lpthread
LIBS += -lz


HEADERS += main.h \
    ArchiveFile.h \
    CompressorModule.h \
    CompressorOperator.h \
    CompressedBlockData.h \
    Params.h \
    FastqCompressor.h \
    ContigBuilder.h \
    ../fastore_bin/BinFile.h \
    ../fastore_rebin/BinFileExtractor.h \
    ../core/Globals.h \
    ../core/FileStream.h \
    ../core/FastqParser.h \
    ../core/FastqPacker.h \
    ../core/FastqRecord.h \
    ../core/DataStream.h \
    ../core/Buffer.h \
    ../core/BitMemory.h \
    ../core/FastqCategorizer.h \
    ../core/Node.h \
    ../core/NodesPacker.h \
    ../core/ReadsClassifier.h \
    ../core/version.h \
    ../ppmd/PPMd.h \
    ../ppmd/Coder.hpp \
    ../ppmd/PPMd.h \
    ../ppmd/PPMdType.h \
    ../ppmd/Stream.hpp \
    ../ppmd/SubAlloc.hpp \
    ../rle/RleEncoder.h \
    ../rc/RangeCoder.h \
    ../rc/SymbolCoderRC.h \
    ../rc/ContextEncoder.h \
    ../qvz/QVZ.h \
    ../qvz/utils.h \
    ../qvz/pmf.h \
    ../qvz/well.h \
    ../qvz/distortion.h \
    ../qvz/quantizer.h \
    ../qvz/codebook.h \
    ../qvz/qv_compressor.h

SOURCES += main.cpp \
    FastqCompressor.cpp \
    ArchiveFile.cpp \
    CompressorModule.cpp \
    CompressorOperator.cpp \
    ContigBuilder.cpp \
    ../fastore_rebin/BinFileExtractor.cpp \
    ../fastore_bin/BinFile.cpp \
    ../core/FileStream.cpp \
    ../core/FastqParser.cpp \
    ../core/FastqPacker.cpp \
    ../core/ReadsClassifier.cpp \
    ../core/FastqCategorizer.cpp \
    ../core/NodesPacker.cpp \
    ../core/version.cpp \
    ../ppmd/Model.cpp \
    ../ppmd/PPMd.cpp \
    ../qvz/Stats.cpp \
    ../qvz/QVZ.cpp \
    ../qvz/pmf.cpp \
    ../qvz/well.cpp \
    ../qvz/distortion.cpp \
    ../qvz/quantizer.cpp \
    ../qvz/codebook.cpp \
    ../qvz/qv_compressor.cpp \
    ../qvz/qv_stream.cpp \
    ../qvz/arith.cpp
