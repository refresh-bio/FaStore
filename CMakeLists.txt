cmake_minimum_required(VERSION 3.0)
project(fastore)

set(FASTORE_VERSION_MAJOR 0)
set(FASTORE_VERSION_MINOR 8)
set(FASTORE_VERSION_PATH 1)


# building prequisites
#
find_package(ZLIB REQUIRED)


# shared compiler configuration
#
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -DNDEBUG -flto")
set(CXX_LD_FLAGS "-lz")

if (APPLE)
	set(CMAKE_CXX_COMPILER "clang++")
	set(CXX_LD_FLAGS "${CXX_LD_FLAGS} -lpthread")
else()
  set(CMAKE_CXX_COMPILER "g++")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -static -fwhole-program")
  set(CXX_LD_FLAGS "${CXX_LD_FLAGS} -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
endif()


# directories setup
#
set(FASTORE_SRC_DIR fastore/)
include_directories(${FASTORE_SRC_DIR})

set(FASTORE_BIN_DIR fastore/fastore_bin)
set(FASTORE_REBIN_DIR fastore/fastore_rebin)
set(FASTORE_PACK_DIR fastore/fastore_pack)

set(QVZ_DIR fastore/fastore_pack)
set(PPMD_DIR fastore/ppmd)


# qvz source files
#
set(QVZ_CORE_SRC ${QVZ_DIR}/codebook.cpp
	${QVZ_DIR}/distortion.cpp
	${QVZ_DIR}/pmf.cpp
	${QVZ_DIR}/quantizer.cpp
	${QVZ_DIR}/util.cpp
	${QVZ_DIR}/well.cpp
	${FASTORE_BIN_DIR}/QVZ.cpp
	)

set(QVZ_COMP_SRC ${QVZ_DIR}/qv_compressor.cpp
	${QVZ_DIR}/arith.cpp 
	${QVZ_DIR}/qv_stream.cpp
	)


# ppmd source files
#
set(PPMD_SRC ${PPMD_DIR}/PPMd.cpp ${PPMD_DIR}/Model.cpp)


# version file setup
#
set(VERSION_FILE version.cpp)
set(VERSION_SRC "${PROJECT_BINARY_DIR}/${VERSION_FILE}")

configure_file( 
	${FASTORE_SRC_DIR}/${VERSION_FILE}.in
	${PROJECT_BINARY_DIR}/${VERSION_FILE}
	)


# fastore_bin configuration
#
#file(GLOB FASTORE_BIN_SRC "fastore/fastore_bin/*.cpp")
set(FASTORE_CORE_SRC ${FASTORE_BIN_DIR}/BinFile.cpp
	${FASTORE_BIN_DIR}/FastqPacker.cpp
	${FASTORE_BIN_DIR}/FastqCategorizer.cpp 
	${FASTORE_BIN_DIR}/FastqParser.cpp
	${FASTORE_BIN_DIR}/FastqStream.cpp
	${FASTORE_BIN_DIR}/FileStream.cpp
	${FASTORE_BIN_DIR}/Stats.cpp
	)

set(FASTORE_BIN_SRC ${FASTORE_BIN_DIR}/main.cpp
	${FASTORE_BIN_DIR}/BinModule.cpp
	${FASTORE_BIN_DIR}/BinOperator.cpp
	${FASTORE_CORE_SRC}
	)

add_executable(fastore_bin ${FASTORE_BIN_SRC} 
	${QVZ_CORE_SRC} 
	${VERSION_SRC})

target_link_libraries(fastore_bin "${CXX_LD_FLAGS}")


# fastore_rebin configuration
#
set(FASTORE_REBIN_SRC ${FASTORE_REBIN_DIR}/main.cpp
	${FASTORE_REBIN_DIR}/NodesPacker.cpp
	${FASTORE_REBIN_DIR}/RebinOperator.cpp
	${FASTORE_REBIN_DIR}/RebinModule.cpp
	${FASTORE_REBIN_DIR}/DnaRebalancer.cpp
	${FASTORE_PACK_DIR}/BinFileExtractor.cpp
	${FASTORE_PACK_DIR}/ContigBuilder.cpp
	${FASTORE_PACK_DIR}/ReadsClassifier.cpp
	${FASTORE_BIN_DIR}/BinOperator.cpp
	${FASTORE_CORE_SRC}
	)

add_executable(fastore_rebin ${FASTORE_REBIN_SRC} 
	${QVZ_CORE_SRC} 
	${VERSION_SRC})

target_link_libraries(fastore_rebin "${CXX_LD_FLAGS}")


# fastore_pack configuration
#
set(FASTORE_PACK_SRC ${FASTORE_PACK_DIR}/main.cpp
	${FASTORE_PACK_DIR}/CompressorModule.cpp
	${FASTORE_PACK_DIR}/CompressorOperator.cpp
	${FASTORE_PACK_DIR}/ArchiveFile.cpp
	${FASTORE_PACK_DIR}/BinFileExtractor.cpp
	${FASTORE_PACK_DIR}/FastqCompressor.cpp
	${FASTORE_PACK_DIR}/ReadsClassifier.cpp
	${FASTORE_PACK_DIR}/ContigBuilder.cpp
	${FASTORE_REBIN_DIR}/NodesPacker.cpp
	${FASTORE_CORE_SRC}
	)

add_executable(fastore_pack ${FASTORE_PACK_SRC} 
	${QVZ_CORE_SRC} ${QVZ_COMP_SRC} 
	${PPMD_SRC} 
	${VERSION_SRC})

target_link_libraries(fastore_pack "${CXX_LD_FLAGS}")
