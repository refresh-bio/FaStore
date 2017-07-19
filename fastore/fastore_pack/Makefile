.PHONY: fastore_pack

all: fastore_pack

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

CXX_OBJS =  CompressorModule.o \
	CompressorOperator.o \
	ArchiveFile.o \
	BinFileExtractor.o \
	FastqCompressor.o \
	ReadsClassifier.o \
	ContigBuilder.o \
	../fastore_bin/BinFile.o \
	../fastore_bin/FastqPacker.o \
	../fastore_bin/FastqParser.o \
	../fastore_bin/FileStream.o \
	../fastore_bin/Stats.o \
	../fastore_bin/FastqCategorizer.o \
	../fastore_rebin/NodesPacker.o

QVZ_OBJS2 = qv_compressor.o arith.o qv_stream.o

PPMD_OBJS = ../ppmd/PPMd.o \
	../ppmd/Model.o

CXX_LIBS += -lz

.cpp.o:
	$(CXX) $(CXX_FLAGS) $(DBG_FLAGS) $(OPT_FLAGS) -c $< -o $@

fastore_pack: $(CXX_OBJS) $(QVZ_OBJS) $(QVZ_OBJS2) $(PPMD_OBJS) mk_version
	$(CXX) $(CXX_FLAGS) $(DBG_FLAGS) $(OPT_FLAGS) $(LD_FLAGS) -o $@ main.cpp version.cpp $(CXX_OBJS) $(QVZ_OBJS) $(QVZ_OBJS2)  $(RLE_OBJS) $(PPMD_OBJS) $(CXX_LIBS)
#	strip $@

mk_version:
	echo "#include \"../fastore_bin/version.h\"" > version.cpp
	echo "std::string GetCompilationTime() {return std::string(__DATE__) + \" -- \" + std::string(__TIME__);}" >> version.cpp
	echo "std::string GetAppVersion() {return \"0.8.0\";}" >> version.cpp

clean:
	-rm -f $(CXX_OBJS)
	-rm -f $(PPMD_OBJS)
	-rm -f $(QVZ_OBJS2)
	-rm -f fastore_pack
