cmake_minimum_required(VERSION 3.16)
project(lact_sim_user LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(include)

find_package(Git QUIET)
set(LACT_PRODUCER_VERSION "source-tree")
if(GIT_FOUND)
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" describe --always --dirty
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE LACT_GIT_DESCRIBE
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(NOT LACT_GIT_DESCRIBE STREQUAL "")
        set(LACT_PRODUCER_VERSION "${LACT_GIT_DESCRIBE}")
    endif()
endif()
message(STATUS "LACT_sim producer version: ${LACT_PRODUCER_VERSION}")

option(LACT_ENABLE_HESSIO "Enable CORSIKA/EventIO support through hessioxxx" ON)
set(HESSIO_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/external/hessioxxx/source" CACHE PATH "hessioxxx source/build root")
if(LACT_ENABLE_HESSIO)
    find_library(HESSIO_LIBRARY hessio PATHS "${HESSIO_ROOT}/lib" NO_DEFAULT_PATH)
else()
    set(HESSIO_LIBRARY "")
endif()
if(HESSIO_LIBRARY)
    set(LACT_HAS_HESSIO TRUE)
    list(APPEND CMAKE_BUILD_RPATH "${HESSIO_ROOT}/lib")
    list(APPEND CMAKE_INSTALL_RPATH "${HESSIO_ROOT}/lib")
else()
    set(LACT_HAS_HESSIO FALSE)
    if(LACT_ENABLE_HESSIO)
        message(WARNING "libhessio not found under ${HESSIO_ROOT}/lib; source.mode=EventIO will be disabled")
    endif()
endif()

option(LACT_ENABLE_HDF5 "Enable native HDF5 output" ON)
if(LACT_ENABLE_HDF5)
    find_package(HDF5 COMPONENTS C)
endif()
if(HDF5_FOUND)
    set(LACT_HAS_HDF5 TRUE)
    include_directories(BEFORE ${HDF5_INCLUDE_DIRS})
else()
    set(LACT_HAS_HDF5 FALSE)
    if(LACT_ENABLE_HDF5)
        message(WARNING "HDF5 C library not found; HDF5 output will be disabled")
    endif()
endif()

option(LACT_ENABLE_ROOT "Enable LACT ROOT output for pyLAST workflows" ON)
if(LACT_ENABLE_ROOT)
    find_package(ROOT 6.24 QUIET COMPONENTS Core RIO Tree)
endif()
if(ROOT_FOUND)
    set(LACT_HAS_ROOT TRUE)
    message(STATUS "Found ROOT ${ROOT_VERSION}; LACT ROOT output enabled")
else()
    set(LACT_HAS_ROOT FALSE)
    if(LACT_ENABLE_ROOT)
        message(STATUS "ROOT >= 6.24 not found; LACT ROOT output disabled")
    endif()
endif()

add_library(lact_sim_core
    src/io/SyntheticPhotonSource.cpp
    src/io/PhotonCsvSource.cpp
    src/io/MirrorFacetCsvReader.cpp
    src/io/EventIOPhotonSemantics.cpp
    src/io/EventIOArrayTiming.cpp
    src/io/Sha256.cpp
    src/app/OpticalSimCommon.cpp
    src/app/ArrayTimingCorrection.cpp
    src/app/PhotonResponseSampler.cpp
    src/app/TelescopeOpticsCache.cpp
    src/app/TriggerResponse.cpp
    src/geometry/FacetFactory.cpp
    src/optics/AtmosphereTransmission.cpp
    src/optics/OpticalTracer.cpp
)
target_compile_definitions(lact_sim_core PUBLIC LACT_PRODUCER_VERSION="${LACT_PRODUCER_VERSION}")

if(LACT_HAS_HESSIO)
    target_sources(lact_sim_core PRIVATE src/io/EventIOPhotonSource.cpp)
    target_include_directories(lact_sim_core PRIVATE "${HESSIO_ROOT}/include")
    target_link_libraries(lact_sim_core PUBLIC "${HESSIO_LIBRARY}")
    target_compile_definitions(lact_sim_core PUBLIC LACT_HAS_HESSIO=1)
endif()

add_executable(run_optical_sim apps/run_optical_sim.cpp)
target_link_libraries(run_optical_sim PRIVATE lact_sim_core)

add_executable(compute_nsb_rate apps/compute_nsb_rate.cpp)
target_link_libraries(compute_nsb_rate PRIVATE lact_sim_core)

if(LACT_HAS_HESSIO)
    add_executable(run_corsika_trace apps/run_corsika_trace.cpp)
    set_target_properties(run_corsika_trace PROPERTIES NO_SYSTEM_FROM_IMPORTED TRUE)
    target_link_libraries(run_corsika_trace PRIVATE lact_sim_core)
    if(LACT_HAS_HDF5)
        target_sources(run_corsika_trace PRIVATE src/io/Hdf5WaveformWriter.cpp)
        target_compile_definitions(run_corsika_trace PRIVATE LACT_HAS_HDF5=1)
        target_link_libraries(run_corsika_trace PRIVATE ${HDF5_LIBRARIES})
    endif()
    if(LACT_HAS_ROOT)
        target_sources(run_corsika_trace PRIVATE src/io/LactEventRootWriter.cpp)
        target_compile_definitions(run_corsika_trace PRIVATE LACT_HAS_ROOT=1)
        target_link_libraries(run_corsika_trace PRIVATE ROOT::Core ROOT::RIO ROOT::Tree)
    endif()
endif()
