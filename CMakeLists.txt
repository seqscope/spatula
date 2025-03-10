cmake_minimum_required(VERSION 2.8...3.23)
project(spatula)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread -std=c++11")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

find_path(HTS_INCLUDE_DIRS htslib/bgzf.h ../htslib/bgzf.h HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/htslib ../htslib)
if(NOT HTS_INCLUDE_DIRS )
    message(FATAL_ERROR "libhts HTS_INCLUDE_DIRS not found")
endif()

message(STATUS ${HTS_INCLUDE_DIRS})

find_library(HTS_LIBRARIES libhts.a HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/htslib/ ../htslib/)
if(NOT HTS_LIBRARIES)
    message(FATAL_ERROR "libhts HTS_LIBRARIES not found")
endif()

message(STATUS ${HTS_LIBRARIES})

find_path(QGEN_INCLUDE_DIRS qgenlib/qgen_error.h  HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/qgenlib ../qgenlib)
if(NOT QGEN_INCLUDE_DIRS )
    message(FATAL_ERROR "libqgen QGEN_INCLUDE_DIRS not found")
endif()

message(STATUS ${QGEN_INCLUDE_DIRS})

find_library(QGEN_LIBRARIES libqgen.a  HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/qgenlib/lib ../qgenlib/lib/)
if(NOT QGEN_LIBRARIES)
    message(FATAL_ERROR "libqgen QGEN_LIBRARIES not found")
endif()

message(STATUS ${QGEN_LIBRARIES})

include_directories(${HTS_INCLUDE_DIRS} ${QGEN_INCLUDE_DIRS})

set(SOURCE_FILES
    cmd_join_pixel_tsv.cpp
    cmd_custom_demux_fastq.cpp
    cmd_stratify_fastq_by_barcodes.cpp
    cmd_filter_common_barcodes.cpp
    cmd_hist.cpp
    cmd_draw_xy.cpp
    cmd_draw_3way.cpp
    cmd_draw_sge.cpp
    cmd_eval_dups_sbcd.cpp
    cmd_combine_sge.cpp
    cmd_combine_sbcd.cpp
    cmd_subset_sge.cpp
    cmd_make_spatial_dge.cpp
    cmd_build_spatial_barcodes.cpp
    cmd_match_tag.cpp
    cmd_search_tag.cpp    
    cmd_match_spatial_barcodes.cpp
    cmd_merge_matched_tags.cpp
    cmd_filter_fastqs.cpp
    cmd_reformat_fastqs.cpp
    cmd_reformat_fastqs_mt.cpp
    cmd_dge2sdge.cpp
    cmd_sge2tsv.cpp
    cmd_convert_sge.cpp
    cmd_segment_sge.cpp
    cmd_filter_tsv.cpp
    multiproc_compressor.h
    sge.cpp
    sge.h
    sge2.cpp
    sge2.h
    tiles.cpp
    tiles.h
    fastq_utils.cpp
    fastq_utils.h
    seq_utils.cpp
    seq_utils.h
    file_utils.cpp
    file_utils.h
    spatula.cpp
    spatula.h
    polygon.h
    polygon.cpp
    nlohmann/json.hpp
    cimg/CImg.h
    )

SET(APP_EXE spatula)

#find_package(OpenMP REQUIRED)
ADD_EXECUTABLE(${APP_EXE} ${SOURCE_FILES})

find_library(ZLIB z HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT ZLIB)
    message(FATAL_ERROR "libz library not found")
endif()

find_library(CURLLIB curl HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CURLLIB)
    message(FATAL_ERROR "libcurl library not found")
endif()

find_package (BZip2)
if (NOT BZIP2_FOUND)
    message(FATAL_ERROR "libbz2 library not found")
else()
    include_directories(${BZIP2_INCLUDE_DIRS})
    set(BZIP2 ${BZIP2_LIBRARIES})
endif()

find_library(LZMA lzma HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT LZMA)
    message(WARNING "liblzma library not found, if you specified --disable-lzma when compiling libhts, please ignore this warning")
    set(LZMA "")
else()
    set(LZMA ${LZMA_LIBRARIES})
endif()

#find_library(DEFLATELIB deflate libdeflate.a HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ../libdeflate/)
#if(NOT DEFLATELIB)
#    message(FATAL_ERROR "libdeflate library not found")
#endif()

find_library(CRYPTOLIB crypto HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CRYPTOLIB)
    message(FATAL_ERROR "libcrypto library not found")
endif()

#TARGET_LINK_LIBRARIES(${APP_EXE} PRIVATE ${QGEN_LIBRARIES} ${HTS_LIBRARIES} ${ZLIB} ${LZMA} ${BZIP2} ${CURLLIB} ${DEFLATELIB} ${CRYPTOLIB} "${OpenMP_CXX_FLAGS}")
#target_compile_options(${APP_EXE} PRIVATE "${OpenMP_CXX_FLAGS}")
#TARGET_LINK_LIBRARIES(${APP_EXE} PRIVATE ${QGEN_LIBRARIES} ${HTS_LIBRARIES} ${ZLIB} ${LZMA} ${BZIP2} ${CURLLIB} ${DEFLATELIB} ${CRYPTOLIB})
TARGET_LINK_LIBRARIES(${APP_EXE} PRIVATE ${QGEN_LIBRARIES} ${HTS_LIBRARIES} ${ZLIB} ${LZMA} ${BZIP2} ${CURLLIB} ${CRYPTOLIB})
target_compile_options(${APP_EXE} PRIVATE)

install(TARGETS ${APP_EXE} RUNTIME DESTINATION bin)
