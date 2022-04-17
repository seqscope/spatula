cmake_minimum_required(VERSION 2.8...3.23)
project(spatula)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread -std=c++11")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

find_path(HTS_INCLUDE_DIRS htslib/bgzf.h ../htslib/bgzf.h HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ../htslib)
if(NOT HTS_INCLUDE_DIRS )
    message(FATAL_ERROR "libhts HTS_INCLUDE_DIRS not found")
endif()

message(STATUS ${HTS_INCLUDE_DIRS})

find_library(HTS_LIBRARIES libhts.a HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ../htslib/)
if(NOT HTS_LIBRARIES)
    message(FATAL_ERROR "libhts HTS_LIBRARIES not found")
endif()

message(STATUS ${HTS_LIBRARIES})

find_path(QGEN_INCLUDE_DIRS qgenlib/qgen_error.h  HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ../qgenlib)
if(NOT QGEN_INCLUDE_DIRS )
    message(FATAL_ERROR "libqgen QGEN_INCLUDE_DIRS not found")
endif()

message(STATUS ${QGEN_INCLUDE_DIRS})

find_library(QGEN_LIBRARIES libqgen.a  HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ../qgenlib/lib/)
if(NOT QGEN_LIBRARIES)
    message(FATAL_ERROR "libqgen QGEN_LIBRARIES not found")
endif()

message(STATUS ${QGEN_LIBRARIES})

include_directories(${HTS_INCLUDE_DIRS} ${QGEN_INCLUDE_DIRS})

set(SOURCE_FILES
    cmd_make_spatial_dge.cpp
    cmd_build_spatial_barcodes.cpp
    cmd_match_spatial_barcodes.cpp    
    seq_utils.cpp
    seq_utils.h
    spatula.cpp
    spatula.h
    )

SET(APP_EXE spatula)

ADD_EXECUTABLE(${APP_EXE} ${SOURCE_FILES}) 

#TARGET_LINK_LIBRARIES(${APP_EXE} ${HTS_LIBRARIES} ${QGEN_LIBRARIES} ${ZLIB} -lhts)
TARGET_LINK_LIBRARIES(${APP_EXE} ${QGEN_LIBRARIES} ${HTS_LIBRARIES} ${ZLIB} -lcurl -llzma -lssl -lcrypto -lbz2)

find_library(ZLIB z HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT ZLIB)
    message(FATAL_ERROR "libz library not found")
endif()

install(TARGETS ${APP_EXE} RUNTIME DESTINATION bin)
