TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
#QMAKE_CXXFLAGS += -DUSE_BOOST_THREAD
QMAKE_CXXFLAGS += -std=c++11 -pthread
CONFIG += c++11

LIBS += -lz
#LIBS += -lboost_thread -lboost_system
LIBS += -lpthread


SOURCES += main.cpp \
    ../fastore_bin/FileStream.cpp \
    ../fastore_bin/FastqStream.cpp \
    ../fastore_bin/BinFile.cpp \
    ../fastore_bin/FastqCategorizer.cpp \
    ../fastore_bin/FastqPacker.cpp \
    ../fastore_bin/BinOperator.cpp \
    ../fastore_bin/FastqParser.cpp \
    ../fastore_pack/BinFileExtractor.cpp \
    ../fastore_pack/ReadsClassifier.cpp \
    ../fastore_bin/version.cpp \
    ../fastore_bin/Stats.cpp \
    RebinOperator.cpp \
    RebinModule.cpp \
    DnaRebalancer.cpp \
    NodesPacker.cpp \
    ../fastore_bin/QVZ.cpp \
    ../fastore_pack/pmf.cpp \
    ../fastore_pack/well.cpp \
    ../fastore_pack/distortion.cpp \
    ../fastore_pack/quantizer.cpp \
    ../fastore_pack/codebook.cpp

HEADERS += \
    ../fastore_bin/FileStream.h \
    ../fastore_bin/DataStream.h \
    ../fastore_bin/Globals.h \
    ../fastore_bin/Buffer.h \
    ../fastore_bin/Utils.h \
    ../fastore_bin/BitMemory.h \
    ../fastore_bin/BinFile.h \
    ../fastore_bin/DataQueue.h \
    ../fastore_bin/DataPool.h \
    ../fastore_bin/Exception.h \
    ../fastore_bin/BinBlockData.h \
    ../fastore_bin/Params.h \
    ../fastore_bin/Thread.h \
    ../fastore_bin/FastqRecord.h \
    ../fastore_bin/FastqCategorizer.h \
    ../fastore_bin/FastqPacker.h \
    ../fastore_bin/BinOperator.h \
    ../fastore_bin/FastqParser.h \
    ../fastore_pack/BinFileExtractor.h \
    ../fastore_pack/ReadsClassifier.h \
    ../fastore_pack/Params.h \
    ../fastore_pack/version.h \
    main.h \
    Params.h \
    RebinModule.h \
    RebinOperator.h \
    DnaRebalancer.h \
    NodesPacker.h \
    ../fastore_bin/QVZ.h \
    ../fastore_pack/pmf.h \
    ../fastore_pack/well.h \
    ../fastore_pack/distortion.h \
    ../fastore_pack/quantizer.h \
    ../fastore_pack/codebook.h

