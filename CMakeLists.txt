cmake_minimum_required(VERSION 2.8...3.23)
project(spatula)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread -std=c++11")

if(APPLE)
    # On macOS, CMake typically finds Apple Clang by default.
    # This block can be used for Apple-specific settings or checks.

    # It's good practice to check if the identified C++ compiler is Clang.
    # CMAKE_CXX_COMPILER_ID is set by CMake after it identifies the compiler.
    # Common IDs for Clang on macOS are "Clang" or "AppleClang".
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
        message(STATUS "Confirmed Clang/AppleClang C++ compiler for macOS: ${CMAKE_CXX_COMPILER}")
        # You can add Clang-specific compiler flags here if necessary, for example:
        # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++") 
        # Note: -stdlib=libc++ is often the default for Clang on macOS.
    else()
        message(WARNING "Expected Clang or AppleClang compiler on macOS, but found C++ compiler ID: '${CMAKE_CXX_COMPILER_ID}' with path '${CMAKE_CXX_COMPILER}'. The project is optimized for Clang; please ensure your environment is configured to use Clang if you encounter issues. If you are using homebrew for example, you may want to run cmake with additional option like '-DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++'. ")
    endif()

    # You can also set macOS specific deployment target if needed:
    # set(CMAKE_OSX_DEPLOYMENT_TARGET "10.13" CACHE STRING "Minimum macOS deployment version")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1 -mpclmul -fno-strict-aliasing")
endif()

# Option to enable finding and linking with libdeflate
option(USE_LIBDEFLATE "Enable if htslib was built with libdeflate support and linking against it is required" OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

find_path(HTS_INCLUDE_DIRS htslib/bgzf.h ../htslib/bgzf.h HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ submodules/htslib ../htslib)
if(NOT HTS_INCLUDE_DIRS )
    message(FATAL_ERROR "libhts HTS_INCLUDE_DIRS not found")
endif()

message(STATUS ${HTS_INCLUDE_DIRS})

find_library(HTS_LIBRARIES hts
  HINTS
    submodules/htslib/
    ../htslib/
    /usr/lib/x86_64-linux-gnu/
    /usr/lib/
    /usr/lib64/
)
if(NOT HTS_LIBRARIES)
    message(FATAL_ERROR "libhts HTS_LIBRARIES not found")
endif()
message(STATUS "HTS: ${HTS_LIBRARIES}")

find_path(QGEN_INCLUDE_DIRS qgenlib/qgen_error.h
  HINTS
    submodules/qgenlib
    ../qgenlib
    /usr/lib/x86_64-linux-gnu/
    /usr/lib/
    /usr/lib64/
)
if(NOT QGEN_INCLUDE_DIRS )
    message(FATAL_ERROR "libqgen QGEN_INCLUDE_DIRS not found")
endif()
message(STATUS "QGEN include: ${QGEN_INCLUDE_DIRS}")

find_library(QGEN_LIBRARIES qgen
  HINTS
    submodules/qgenlib/lib
    ../qgenlib/lib/
    /usr/lib/x86_64-linux-gnu/
    /usr/lib/
    /usr/lib64/
)

if(NOT QGEN_LIBRARIES)
    message(FATAL_ERROR "libqgen QGEN_LIBRARIES not found")
endif()
message(STATUS "QGEN lib: ${QGEN_LIBRARIES}")

include_directories(${HTS_INCLUDE_DIRS} ${QGEN_INCLUDE_DIRS})

set(SOURCE_FILES
    cmd_mex_subset.cpp
    cmd_sptsv2mex.cpp
    cmd_sptsv2model.cpp
    cmd_mex2sptsv.cpp
    cmd_png_mono2rgba.cpp
    cmd_join_pixel_decode.cpp
    cmd_append_topk_tsv.cpp
    cmd_reformat_tsv.cpp
    cmd_pseudobulk_from_decode.cpp
    cmd_diffexp_decode.cpp
    cmd_diffexp_model_matrix.cpp
    cmd_merge_pseudobulk.cpp
    cmd_paste_pixel_tsv.cpp
    cmd_join_pixel_tsv.cpp
    cmd_custom_demux_fastq.cpp
    cmd_stratify_fastq_by_barcodes.cpp
    cmd_filter_common_barcodes.cpp
    cmd_hist.cpp
    cmd_draw_xy.cpp
    cmd_draw_3way.cpp
    cmd_draw_sge.cpp
    cmd_draw_tsv.cpp
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
    generic_utils.cpp
    generic_utils.h
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
    nanoflann/nanoflann.hpp
    fpng/fpng.h
    fpng/fpng.cpp
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

set(LIBDEFLATE_LIBRARY "") # Initialize variable

if(USE_LIBDEFLATE)
    find_library(DEFLATELIB_PATH deflate libdeflate.a HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/ ../libdeflate/)
    if(NOT DEFLATELIB_PATH)
        message(WARNING "USE_LIBDEFLATE is ON, but libdeflate (libdeflate.a, libdeflate.so, libdeflate.dylib, etc.) was not found. "
        "If your htslib requires it, linking may fail or you might encounter runtime errors. "
        "Ensure libdeflate is installed in a standard location, or provide its path via CMAKE_PREFIX_PATH or LIBDEFLATE_DIR (environment variable pointing to the installation prefix).")
    else()
        set(DEFLATELIB ${DEFLATELIB_PATH})
        message(STATUS "Found libdeflate for linking: ${DEFLATELIB}")
    endif()
endif()


find_library(CRYPTOLIB crypto HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
if(NOT CRYPTOLIB)
    message(FATAL_ERROR "libcrypto library not found")
endif()

# # add PNG library
# find_library(PNGLIB png HINTS /usr/lib/x86_64-linux-gnu/ /usr/lib/ /usr/lib64/)
# if(NOT PNGLIB)
#     message(FATAL_ERROR "libpng library not found")
# endif()

#TARGET_LINK_LIBRARIES(${APP_EXE} PRIVATE ${QGEN_LIBRARIES} ${HTS_LIBRARIES} ${ZLIB} ${LZMA} ${BZIP2} ${CURLLIB} ${DEFLATELIB} ${CRYPTOLIB} "${OpenMP_CXX_FLAGS}")
#target_compile_options(${APP_EXE} PRIVATE "${OpenMP_CXX_FLAGS}")
#TARGET_LINK_LIBRARIES(${APP_EXE} PRIVATE ${QGEN_LIBRARIES} ${HTS_LIBRARIES} ${ZLIB} ${LZMA} ${BZIP2} ${CURLLIB} ${CRYPTOLIB})
TARGET_LINK_LIBRARIES(${APP_EXE} PRIVATE ${QGEN_LIBRARIES} ${HTS_LIBRARIES} ${ZLIB} ${LZMA} ${BZIP2} ${CURLLIB} ${CRYPTOLIB} ${DEFLATELIB})
target_compile_options(${APP_EXE} PRIVATE)

install(TARGETS ${APP_EXE} RUNTIME DESTINATION bin)
