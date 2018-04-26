TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

QMAKE_CXXFLAGS += -std=c++11
QMAKE_CXXFLAGS += -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE

LIBS += -lz
LIBS += -lpthread


HEADERS += main.h \
    BinModule.h \
    BinOperator.h \
    BinFile.h \
    Params.h \
    ../core/BinBlockData.h \
    ../core/FileStream.h \
    ../core/FastqStream.h \
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
    ../core/Node.h \
    ../core/version.h \
    ../qvz/Quality.h \
    ../qvz/Stats.h \
    ../qvz/QVZ.h \
    ../qvz/pmf.h \
    ../qvz/well.h \
    ../qvz/distortion.h \
    ../qvz/quantizer.h \
    ../qvz/codebook.h

SOURCES += main.cpp \
    BinModule.cpp \
    BinOperator.cpp \
    BinFile.cpp \
    ../core/FileStream.cpp \
    ../core/FastqStream.cpp \
    ../core/FastqCategorizer.cpp \
    ../core/FastqPacker.cpp \
    ../core/FastqParser.cpp \
    ../core/version.cpp \
    ../qvz/pmf.cpp \
    ../qvz/well.cpp \
    ../qvz/distortion.cpp \
    ../qvz/quantizer.cpp \
    ../qvz/codebook.cpp \
    ../qvz/QVZ.cpp \
    ../qvz/Stats.cpp
