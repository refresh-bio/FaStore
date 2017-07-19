TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
#QMAKE_CXXFLAGS += -DUSE_BOOST_THREAD

CONFIG += c++11
QMAKE_CXXFLAGS += -std=c++11

LIBS += -lpthread
LIBS += -lz
#LIBS += -lboost_thread -lboost_system


HEADERS += \
    ../fastore_bin/utils.h \
    ../fastore_bin/Globals.h \
    ../fastore_bin/FileStream.h \
    ../fastore_bin/FastqParser.h \
    ../fastore_bin/FastqPacker.h \
    ../fastore_bin/FastqRecord.h \
    ../fastore_bin/DataStream.h \
    ../fastore_bin/Buffer.h \
    ../fastore_bin/BitMemory.h \
    ../fastore_bin/BinFile.h \
    ../fastore_bin/FastqCategorizer.h \
    ../fastore_bin/version.h \
    ../fastore_bin/Node.h \
    ../fastore_rebin/NodesPacker.h \
    BinFileExtractor.h \
    ArchiveFile.h \
    CompressorModule.h \
    CompressorOperator.h \
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
    CompressedBlockData.h \
    Params.h \
    main.h \
    FastqCompressor.h \
    ReadsClassifier.h \
    ContigBuilder.h \
    ../fastore_bin/QVZ.h \
    pmf.h \
    well.h \
    distortion.h \
    quantizer.h \
    codebook.h \
    qv_compressor.h

SOURCES += \
    main.cpp \
    ../fastore_bin/FileStream.cpp \
    ../fastore_bin/FastqParser.cpp \
    ../fastore_bin/FastqPacker.cpp \
    ../fastore_bin/BinFile.cpp \
    ../fastore_bin/Stats.cpp \
    ../fastore_bin/FastqCategorizer.cpp \
    ../fastore_rebin/NodesPacker.cpp \
    version.cpp \
    BinFileExtractor.cpp \
    FastqCompressor.cpp \
    ArchiveFile.cpp \
    CompressorModule.cpp \
    CompressorOperator.cpp \
    ../ppmd/Model.cpp \
    ../ppmd/PPMd.cpp \
    ReadsClassifier.cpp \
    ContigBuilder.cpp \
    ../fastore_bin/QVZ.cpp \
    pmf.cpp \
    well.cpp \
    distortion.cpp \
    quantizer.cpp \
    codebook.cpp \
    qv_compressor.cpp \
    qv_stream.cpp \
    arith.cpp
