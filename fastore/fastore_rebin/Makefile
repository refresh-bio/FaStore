.PHONY: fastore_bin

all: fastore_rebin

ifndef CXX
CXX = g++
endif

ifndef DBG_FLAGS
DBG_FLAGS += -DNDEBUG
endif

ifndef OPT_FLAGS
OPT_FLAGS += -O2 -flto -fwhole-program
endif

CXX_FLAGS += -m64 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
CXX_FLAGS += -std=c++11 -pthread

#LD_FLAGS += -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
#LD_FLAGS += -static

CXX_OBJS = ../fastore_bin/FileStream.o \
    ../fastore_bin/FastqStream.o \
    ../fastore_bin/BinFile.o \
    ../fastore_bin/FastqCategorizer.o \
    ../fastore_bin/FastqPacker.o \
    ../fastore_bin/BinOperator.o \
    ../fastore_bin/FastqParser.o \
    ../fastore_bin/Stats.o \
    ../fastore_bin/version.o \
    ../fastore_pack/BinFileExtractor.o \
    ../fastore_pack/ContigBuilder.o \
    ../fastore_pack/ReadsClassifier.o \
    NodesPacker.o \
    RebinOperator.o \
    RebinModule.o \
    DnaRebalancer.o


CXX_LIBS += -lz

.cpp.o:
	$(CXX) $(CXX_FLAGS) $(DBG_FLAGS) $(OPT_FLAGS) -c $< -o $@

fastore_rebin: $(CXX_OBJS) $(QVZ_OBJS)
	$(CXX) $(CXX_FLAGS) $(DBG_FLAGS) $(OPT_FLAGS) $(LD_FLAGS) -o $@ main.cpp $(CXX_OBJS) $(QVZ_OBJS) $(CXX_LIBS)
#	strip $@

clean:
	-rm -f $(CXX_OBJS)
	-rm -f fastore_rebin
