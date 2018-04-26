TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

QMAKE_CXXFLAGS += -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
QMAKE_CXXFLAGS += -std=c++11 -pthread

LIBS += -lz
LIBS += -lpthread


HEADERS += main.h \
    Params.h \
    RebinModule.h \
    RebinOperator.h \
    DnaRebalancer.h \
    BinFileExtractor.h \
    ../fastore_bin/BinFile.h \
    ../fastore_bin/BinOperator.h \
    ../fastore_bin/Params.h \
    ../fastore_pack/Params.h \
    ../core/BinBlockData.h \
    ../core/NodesPacker.h \
    ../core/FileStream.h \
    ../core/DataStream.h \
    ../core/Globals.h \
    ../core/Buffer.h \
    ../core/Utils.h \
    ../core/BitMemory.h \
    ../core/DataQueue.h \
    ../core/DataPool.h \
    ../core/Exception.h \
    ../core/Thread.h \
    ../core/FastqRecord.h \
    ../core/FastqCategorizer.h \
    ../core/FastqPacker.h \
    ../core/FastqParser.h \
    ../core/ReadsClassifier.h \
    ../core/version.h \
    ../qvz/QVZ.h \
    ../qvz/pmf.h \
    ../qvz/well.h \
    ../qvz/distortion.h \
    ../qvz/quantizer.h \
    ../qvz/codebook.h

SOURCES += main.cpp \
    BinFileExtractor.cpp \
    RebinOperator.cpp \
    RebinModule.cpp \
    DnaRebalancer.cpp \
    ../fastore_bin/BinFile.cpp \
    ../fastore_bin/BinOperator.cpp \
    ../core/FileStream.cpp \
    ../core/FastqStream.cpp \
    ../core/FastqCategorizer.cpp \
    ../core/FastqPacker.cpp \
    ../core/FastqParser.cpp \
    ../core/ReadsClassifier.cpp \
    ../core/NodesPacker.cpp \
    ../core/version.cpp \
    ../qvz/Stats.cpp \
    ../qvz/QVZ.cpp \
    ../qvz/pmf.cpp \
    ../qvz/well.cpp \
    ../qvz/distortion.cpp \
    ../qvz/quantizer.cpp \
    ../qvz/codebook.cpp
