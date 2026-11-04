cmake_minimum_required(VERSION 3.20...3.31)

project(rn_obj_search CXX)

# Select a compiler with -DCMAKE_CXX_COMPILER=... (or the CXX environment variable) when configuring
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
# Note: NDEBUG is deliberately left undefined in release builds. The project's own checks do not depend on it (debug.h), but the
# asserts of the third-party headers do, and defining it would change the hash-table code (sparsehash) that the queries run
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
# Pass -DRN_NATIVE_ARCH=OFF to build a binary that runs on other CPUs
# Note: -ffp-contract=off prevents FMA from changing floating-point results (Clang contracts by default)
option(RN_NATIVE_ARCH "Optimise for the build machine's CPU (binary may not run on other CPUs)" ON)
if (RN_NATIVE_ARCH)
    string(APPEND CMAKE_CXX_FLAGS_RELEASE " -march=native -ffp-contract=off")
endif()

# Dependencies (checked at configure time so a missing package is reported clearly)
find_package(Boost 1.70 CONFIG COMPONENTS serialization)
if (NOT Boost_FOUND)
    message(FATAL_ERROR "Boost.Serialization (>= 1.70) not found. On Ubuntu/Debian install libboost-serialization-dev, "
                        "or pass -DBoost_DIR=/path/to/lib/cmake/Boost-<version>")
endif()

find_path(METIS_INCLUDE_DIR metis.h)
find_library(METIS_LIBRARY metis)
if (NOT METIS_INCLUDE_DIR OR NOT METIS_LIBRARY)
    message(FATAL_ERROR "METIS not found. On Ubuntu/Debian install libmetis-dev, "
                        "or pass -DMETIS_INCLUDE_DIR=... and -DMETIS_LIBRARY=...")
endif()
add_library(METIS::METIS UNKNOWN IMPORTED)
set_target_properties(METIS::METIS PROPERTIES
    IMPORTED_LOCATION "${METIS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIR}")

# Sparsehash is header-only
find_path(SPARSEHASH_INCLUDE_DIR google/dense_hash_map)
if (NOT SPARSEHASH_INCLUDE_DIR)
    message(FATAL_ERROR "Google Sparsehash not found. On Ubuntu/Debian install libsparsehash-dev, "
                        "or pass -DSPARSEHASH_INCLUDE_DIR=...")
endif()

# Needed by the verification's ground truth computation and by rn_obj_search_par (below)
find_package(OpenMP)

# CONFIGURE_DEPENDS re-globs on each build, so new source files are picked up without re-running cmake
file(GLOB MAIN_SRC CONFIGURE_DEPENDS
    "main.cpp"
    "command/*.cpp"
    "tuple/*.cpp"
    "utility/*.cpp"
    "queue/*.cpp"
    "processing/*.cpp"
    "verification/*.cpp"
)

# Standard Optimisations
# By default all optimisations are already applied without pre-processor directives
add_executable(rn_obj_search ${MAIN_SRC})

# Additional Executables for Experiments

# Version that collects other query statistics (e.g. nodes visited in INE), and also runs the debug-only checks
# (RN_DEBUG_ASSERT in debug.h), as its times are never used
# Note: This is not to be used for query times due to collection overhead
add_executable(rn_obj_search_stats ${MAIN_SRC})
target_compile_definitions(rn_obj_search_stats PRIVATE COLLECT_STATISTICS RN_DEBUG_CHECKS)

# Version that compiles the OpenMP pragmas written as PARALLEL_BUILDS_PRAGMA (common.h), used by the parameter-test SUL-Tree
# builds (runParameterExperiments build_rn): several SUL-Trees at once, each processing its tree nodes in parallel
# Note: This is not to be used for query times or reported index construction times
add_executable(rn_obj_search_par ${MAIN_SRC})
target_compile_definitions(rn_obj_search_par PRIVATE PARALLEL_BUILDS)

foreach(target rn_obj_search rn_obj_search_stats rn_obj_search_par)
    # Debug-only checks (RN_DEBUG_ASSERT in debug.h) in the Debug configuration of every executable
    target_compile_definitions(${target} PRIVATE $<$<CONFIG:Debug>:RN_DEBUG_CHECKS>)
    target_include_directories(${target} PRIVATE "${SPARSEHASH_INCLUDE_DIR}")
    target_link_libraries(${target} PRIVATE Boost::serialization METIS::METIS)
    if (OpenMP_CXX_FOUND)
        target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
    endif()
endforeach()
