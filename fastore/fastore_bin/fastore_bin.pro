TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
#QMAKE_CXXFLAGS += -DUSE_BOOST_THREAD

CONFIG += c++11
QMAKE_CXXFLAGS += -std=c++11

LIBS += -lz
#LIBS += -lboost_thread -lboost_system
LIBS += -lpthread


SOURCES += main.cpp \
    FileStream.cpp \
    FastqStream.cpp \
    BinFile.cpp \
    BinModule.cpp \
    BinOperator.cpp \
    FastqCategorizer.cpp \
    FastqPacker.cpp \
    FastqParser.cpp \
    version.cpp \
    ../fastore_pack/pmf.cpp \
    ../fastore_pack/well.cpp \
    ../fastore_pack/distortion.cpp \
    ../fastore_pack/quantizer.cpp \
    ../fastore_pack/codebook.cpp \
    QVZ.cpp \
    Stats.cpp


HEADERS += \
    FileStream.h \
    FastqStream.h \
    DataStream.h \
    Globals.h \
    Buffer.h \
    Utils.h \
    BitMemory.h \
    BinFile.h \
    BinModule.h \
    DataQueue.h \
    DataPool.h \
    BinOperator.h \
    Exception.h \
    BinBlockData.h \
    Params.h \
    Thread.h \
    main.h \
    FastqRecord.h \
    FastqCategorizer.h \
    FastqPacker.h \
    FastqParser.h \
    version.h \
    Node.h \
    Quality.h \
    Stats.h \
    QVZ.h \
    ../fastore_pack/pmf.h \
    ../fastore_pack/well.h \
    ../fastore_pack/distortion.h \
    ../fastore_pack/quantizer.h \
    ../fastore_pack/codebook.h
