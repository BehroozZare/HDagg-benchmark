# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
    set(HDAGG_TOPLEVEL_PROJECT OFF)
else()
    set(HDAGG_TOPLEVEL_PROJECT ON)
endif()


# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.14")
set(CMAKE_POLICY_DEFAULT_CMP0127 NEW)

if(HDAGG_TOPLEVEL_PROJECT)
    cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else()
    # Don't use cmake_minimum_required here to avoid implicitly overriding parent policies
    if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
        message(FATAL_ERROR "CMake required version to build IPC is ${REQUIRED_CMAKE_VERSION}")
    endif()
endif()



list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

project(HDAGG DESCRIPTION "A upgraded version of LBC library with HDagg functionality" LANGUAGES CXX)

################################################################################

#Project options
option(HDAGG_WITH_SPMP "Use SpMP library with HDAGG" ON)
option(HDAGG_WITH_DAGP "Add DAGP for benchmarking" OFF)

option(HDAGG_WITH_DEMO "Add demo executables" ON)



################################################################################

### Configuration
set(HDAGG_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(${PROJECT_NAME}Dependencies)

################################################################################
# LBC Library
################################################################################

# project source files
file(GLOB SRCFILES
        "src/dense_blas/*.cpp"
        "src/hdagg/*.cpp"
        "src/Group/*.cpp"
        "src/lbc/*.cpp"
        "src/profiler/*.cpp"
        "src/sparse_blas/*.cpp"
        "src/utils/*.cpp"
        "src/utils/external/MatrixMarketIO/*.cpp"
        )

add_library(${PROJECT_NAME}_lib ${SRCFILES})

target_include_directories(${PROJECT_NAME}_lib PUBLIC
        "src/dense_blas/includes"
        "src/hdagg/includes"
        "src/Group/includes"
        "src/lbc/includes"
        "src/profiler/includes"
        "src/sparse_blas/includes"
        "src/utils/includes"
        "src/utils/external/MatrixMarketIO/includes"
        )

################################################################################
# Compiler options
################################################################################

add_compile_options(-Wno-write-strings)
add_compile_options(-Wunused-result)


find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME}_lib PUBLIC OpenMP::OpenMP_CXX)

# Figure out SSE level support
message(STATUS "Seaching for SSE...")
find_package(SSE)

# Figure out AVX level support
message(STATUS "Searching for AVX...")
find_package(AVX)

# Figure out FMA level support
message(STATUS "Searching for FMA...")
find_package(FMA)

# Add SSE, AVX, and FMA flags to compiler flags
string(REPLACE " " ";" SIMD_FLAGS "${SSE_FLAGS} ${AVX_FLAGS} ${FMA_FLAGS}")
target_compile_options(${PROJECT_NAME}_lib PUBLIC ${SIMD_FLAGS})

# Add -pthread to compilation and linking
find_package(Threads)
target_link_libraries(${PROJECT_NAME}_lib PUBLIC ${CMAKE_THREAD_LIBS_INIT})

# Use C++11
target_compile_features(${PROJECT_NAME}_lib PUBLIC cxx_std_17)

################################################################################
# Required Libraries
################################################################################

### Find Metis
find_package(METIS OPTIONAL_COMPONENTS)
if(EXISTS "${METIS_INCLUDES}")
    add_definitions(-DMETIS)
    target_link_libraries(${PROJECT_NAME}_lib PUBLIC
            ${METIS_LIBRARY}
            ${METIS_LIBRARIES}
            )
else()
    message(STATUS "The METIS is not found")
endif()

find_package(PAPI OPTIONAL_COMPONENTS)
if(EXISTS "${PAPI_INCLUDE_DIRS}")
    target_include_directories(${PROJECT_NAME}_lib PUBLIC ${PAPI_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME}_lib PUBLIC
            ${PAPI_LIBRARIES}
            )
else()
    message(STATUS "PAPI is not found.")
endif()

### Find Intel MKL if available
find_package(MKL OPTIONAL_COMPONENTS)
if(EXISTS "${MKL_INCLUDE_DIRS}")
    target_include_directories(${PROJECT_NAME}_lib PUBLIC
            ${MKL_INCLUDE_DIRS}
            )
    target_link_libraries(${PROJECT_NAME}_lib PUBLIC
            ${MKL_GNU_THREAD}
            ${PTHREAD}
            ${IOMP5}
            ${MKL_SEQUENTIAL}
            ${MKL_INTEL_THREAD}
            ${MKL_CORE}
            ${M}
            ${DL}
            )

else()
    message(STATUS "Intel MKL not found, building without MKL.")
    set(MKL 0)
endif()


### Check whether DAGP root exists - Relative works
if(HDAGG_WITH_DAGP)
    if(EXISTS "$ENV{DAGPROOT}")
        message(STATUS "The DAGP is found - $ENV{DAGPROOT}.")
        target_include_directories(${PROJECT_NAME}_lib PUBLIC
                "$ENV{DAGPROOT}/src/recBisection/"
                "$ENV{DAGPROOT}/src/common/"
                "src/utils/external/DAGPInterface/includes"
                )
        target_link_libraries(${PROJECT_NAME}_lib PUBLIC
                "$ENV{DAGPROOT}/lib/libdagp.a"
                "src/utils/external/DAGPInterface/*.cpp"
                )
        add_definitions(-DDAGP)
    else()
        message(STATUS "Please defined the DAGP directory using DAGP_ROOT variable")
    endif()
endif()

### Check whether SPMP root exists - Relative works
if(HDAGG_WITH_SPMP)
    if(EXISTS "$ENV{SPMPROOT}")
        message(STATUS "Found SpMP - $ENV{SPMPROOT}.")
        target_include_directories(${PROJECT_NAME}_lib PUBLIC
                $ENV{SPMPROOT}
                "src/utils/external/SpMPInterface/includes"
                )
        target_link_libraries(${PROJECT_NAME}_lib PUBLIC
                "$ENV{SPMPROOT}/libspmp.a"
                )
        add_definitions(-DSPMP)
    else()
        message(STATUS "Please defined the SpMP directory using SPMP_ROOT variable - not found")
    endif()
endif()

# Logger
#target_link_libraries(${PROJECT_NAME}_lib PUBLIC spdlog::spdlog)

# Unit tests
if(HDAGG_WITH_DEMO)
    add_subdirectory(demo)
endif()