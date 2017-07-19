.PHONY: fastore_bin

all: fastore_bin

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

CXX_OBJS = BinModule.o \
	BinOperator.o \
	BinFile.o \
	FastqPacker.o \
	FastqCategorizer.o \
	FastqParser.o \
	FastqStream.o \
	FileStream.o \
	Stats.o

CXX_LIBS += -lz

.cpp.o:
	$(CXX) $(CXX_FLAGS) $(DBG_FLAGS) $(OPT_FLAGS) -c $< -o $@

fastore_bin: $(CXX_OBJS) $(QVZ_OBJS) mk_version	
	$(CXX) $(CXX_FLAGS) $(DBG_FLAGS) $(OPT_FLAGS) $(LD_FLAGS) -o $@ main.cpp version.cpp $(CXX_OBJS) $(QVZ_OBJS) $(CXX_LIBS)
#	strip $@

mk_version:
	echo "#include \"version.h\"" > version.cpp
	echo "std::string GetCompilationTime() {return std::string(__DATE__) + \" -- \" + std::string(__TIME__);}" >> version.cpp
	echo "std::string GetAppVersion() {return \"0.8.0\";}" >> version.cpp

clean:
	-rm -f $(CXX_OBJS)
	-rm -f fastore_bin
