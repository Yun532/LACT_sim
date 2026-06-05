cmake_minimum_required(VERSION 3.16)
project(lact_sim_user LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(include)

option(LACT_ENABLE_HESSIO "Enable direct CORSIKA/EventIO support through hessioxxx" ON)
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
    else()
        message(STATUS "hessio/EventIO support disabled by LACT_ENABLE_HESSIO=OFF")
    endif()
endif()

option(LACT_ENABLE_HDF5 "Enable native HDF5 output for trace programs" ON)
if(LACT_ENABLE_HDF5)
    find_package(HDF5 COMPONENTS C)
endif()
if(HDF5_FOUND)
    set(LACT_HAS_HDF5 TRUE)
else()
    set(LACT_HAS_HDF5 FALSE)
    if(LACT_ENABLE_HDF5)
        message(WARNING "HDF5 C library not found; native .h5 output will be disabled")
    else()
        message(STATUS "native HDF5 output disabled by LACT_ENABLE_HDF5=OFF")
    endif()
endif()

add_library(lact_sim_core
    src/io/SyntheticPhotonSource.cpp
    src/io/PhotonCsvSource.cpp
    src/io/MirrorFacetCsvReader.cpp
    src/app/OpticalSimCommon.cpp
    src/geometry/FacetFactory.cpp
    src/optics/AtmosphereTransmission.cpp
    src/optics/OpticalTracer.cpp
)

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
    target_link_libraries(run_corsika_trace PRIVATE lact_sim_core)
    if(LACT_HAS_HDF5)
        target_compile_definitions(run_corsika_trace PRIVATE LACT_HAS_HDF5=1)
        target_include_directories(run_corsika_trace PRIVATE ${HDF5_INCLUDE_DIRS})
        target_link_libraries(run_corsika_trace PRIVATE ${HDF5_LIBRARIES})
    endif()
endif()
