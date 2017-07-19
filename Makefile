.PHONY: cpp11

all: cpp11

BIN_DIR = bin

export
#DBG_FLAGS += -DDEBUG -DDEV_DEBUG_MODE
DBG_FLAGS += -DNDEBUG
OPT_FLAGS += -O2 
#OPT_FLAGS += -DNDEBUG
#DBG_FLAGS += -g
#OPT_FLAGS += -O3
#OPT_FLAGS += -O0 -g -gdwarf-2 
#DBG_FLAGS += -pg

LD_FLAGS += -static
LD_FLAGS += -Wl,--whole-archive -lpthread -Wl,--no-whole-archive

CXX = g++

QVZ_DIR := $(shell pwd)/fastore/fastore_pack
QVZ_OBJS = $(QVZ_DIR)/codebook.o \
	$(QVZ_DIR)/distortion.o \
	$(QVZ_DIR)/pmf.o \
	$(QVZ_DIR)/quantizer.o \
	$(QVZ_DIR)/util.o \
	$(QVZ_DIR)/well.o \
	$(shell pwd)/fastore/fastore_bin/QVZ.o


cpp11:
	cd fastore/fastore_bin && make fastore_bin
	cd fastore/fastore_rebin && make fastore_rebin
	cd fastore/fastore_pack && make fastore_pack
	test -d $(BIN_DIR) || mkdir $(BIN_DIR)	
	mv fastore/fastore_bin/fastore_bin $(BIN_DIR)/fastore_bin
	mv fastore/fastore_rebin/fastore_rebin $(BIN_DIR)/fastore_rebin
	mv fastore/fastore_pack/fastore_pack $(BIN_DIR)/fastore_pack
	strip $(BIN_DIR)/fastore_*

clean:
	cd fastore/fastore_bin/ && make clean
	cd fastore/fastore_pack/ && make clean
	cd fastore/fastore_rebin/ && make clean
	-rm -rf $(QVZ_OBJS)
	-rm -rf $(BIN_DIR)
