.PHONY: build

CMAKE_PRESENT := $(shell command -v cmake 2> /dev/null)

all: build

BIN_DIR = bin
BUILD_DIR = __build

build:
ifndef CMAKE_PRESENT
    $(error "cmake is not available -- please install cmake")
endif
	test -d $(BUILD_DIR) || mkdir $(BUILD_DIR)	
	cd $(BUILD_DIR) && cmake .. && make
	test -d $(BIN_DIR) || mkdir $(BIN_DIR)	
	cp $(BUILD_DIR)/fastore_* $(BIN_DIR)/
	strip $(BIN_DIR)/fastore_*

clean:
	-rm -rf $(BUILD_DIR)
	-rm -rf $(BIN_DIR)
