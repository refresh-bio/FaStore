all: rle

CXX = g++

CXX_FLAGS += -O3 -DNDEBUG -flto -fwhole-program
CXX_FLAGS += -m64 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -static

rle:
	$(CXX) $(CXX_FLAGS) -o $@ main.cpp rle.cpp
	strip $@

clean:
	-rm rle
