cmake_minimum_required(VERSION 3.16)
project(simlite LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(include)
enable_testing()

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
    include_directories(BEFORE ${HDF5_INCLUDE_DIRS})
else()
    set(LACT_HAS_HDF5 FALSE)
    if(LACT_ENABLE_HDF5)
        message(WARNING "HDF5 C library not found; native .h5 output will be disabled")
    else()
        message(STATUS "native HDF5 output disabled by LACT_ENABLE_HDF5=OFF")
    endif()
endif()

option(LACT_ENABLE_ROOT "Enable ROOT lact_event output for pylast adapter workflows" ON)
if(LACT_ENABLE_ROOT)
    find_package(ROOT 6.24 QUIET COMPONENTS Core RIO Tree)
endif()
if(ROOT_FOUND)
    set(LACT_HAS_ROOT TRUE)
    message(STATUS "Found ROOT ${ROOT_VERSION}; lact_event ROOT output enabled")
else()
    set(LACT_HAS_ROOT FALSE)
    if(LACT_ENABLE_ROOT)
        message(STATUS "ROOT >= 6.24 not found; lact_event ROOT output disabled")
    else()
        message(STATUS "lact_event ROOT output disabled by LACT_ENABLE_ROOT=OFF")
    endif()
endif()

add_library(simlite_core
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
    src/electronics/DetectorPipeline.cpp
    src/geometry/FacetFactory.cpp
    src/optics/AtmosphereTransmission.cpp
    src/optics/OpticalTracer.cpp
)
target_compile_definitions(
    simlite_core PUBLIC
    LACT_PRODUCER_VERSION="${LACT_PRODUCER_VERSION}"
)

if(LACT_HAS_HESSIO)
    target_sources(simlite_core PRIVATE src/io/EventIOPhotonSource.cpp)
    target_include_directories(simlite_core PRIVATE "${HESSIO_ROOT}/include")
    target_link_libraries(simlite_core PUBLIC "${HESSIO_LIBRARY}")
    target_compile_definitions(simlite_core PUBLIC LACT_HAS_HESSIO=1)
endif()

add_executable(test_synthetic_source
    apps/test_synthetic_source.cpp
)
target_link_libraries(test_synthetic_source PRIVATE simlite_core)
add_test(NAME test_synthetic_source COMMAND test_synthetic_source)

add_executable(test_sha256
    apps/test_sha256.cpp
)
target_link_libraries(test_sha256 PRIVATE simlite_core)
add_test(NAME test_sha256 COMMAND test_sha256)

add_executable(run_optical_sim
    apps/run_optical_sim.cpp
)
target_link_libraries(run_optical_sim PRIVATE simlite_core)

add_executable(compute_nsb_rate
    apps/compute_nsb_rate.cpp
)
target_link_libraries(compute_nsb_rate PRIVATE simlite_core)

add_executable(run_camera_electronics
    apps/run_camera_electronics.cpp
)
target_link_libraries(run_camera_electronics PRIVATE simlite_core)

if(LACT_HAS_HESSIO)
    add_executable(run_corsika_trace
        apps/run_corsika_trace.cpp
    )
    set_target_properties(run_corsika_trace PROPERTIES NO_SYSTEM_FROM_IMPORTED TRUE)
    target_link_libraries(run_corsika_trace PRIVATE simlite_core)
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

    add_executable(test_eventio_wavelength
        apps/test_eventio_wavelength.cpp
    )
    target_link_libraries(test_eventio_wavelength PRIVATE simlite_core)
    add_test(NAME test_eventio_wavelength COMMAND test_eventio_wavelength)
endif()

add_executable(test_facet_geometry
    apps/test_facet_geometry.cpp
)
target_link_libraries(test_facet_geometry PRIVATE simlite_core)
add_test(NAME test_facet_geometry COMMAND test_facet_geometry)

add_executable(test_mirror_facet_csv
    apps/test_mirror_facet_csv.cpp
)
target_link_libraries(test_mirror_facet_csv PRIVATE simlite_core)
add_test(NAME test_mirror_facet_csv COMMAND test_mirror_facet_csv)

add_executable(test_optical_surface
    apps/test_optical_surface.cpp
)
target_link_libraries(test_optical_surface PRIVATE simlite_core)
add_test(NAME test_optical_surface COMMAND test_optical_surface)

add_executable(test_raytrace_primitives
    apps/test_raytrace_primitives.cpp
)
target_link_libraries(test_raytrace_primitives PRIVATE simlite_core)
add_test(NAME test_raytrace_primitives COMMAND test_raytrace_primitives)

add_executable(test_dc_center_rays
    apps/test_dc_center_rays.cpp
)
target_link_libraries(test_dc_center_rays PRIVATE simlite_core)
add_test(NAME test_dc_center_rays COMMAND test_dc_center_rays)

add_executable(test_mirror_1229_layout
    apps/test_mirror_1229_layout.cpp
)
target_link_libraries(test_mirror_1229_layout PRIVATE simlite_core)
add_test(NAME test_mirror_1229_layout
         COMMAND test_mirror_1229_layout ${CMAKE_CURRENT_SOURCE_DIR}/configs/mirror_1229_facets.csv)

add_executable(test_camera_geometry
    apps/test_camera_geometry.cpp
)
target_link_libraries(test_camera_geometry PRIVATE simlite_core)
add_test(NAME test_camera_geometry COMMAND test_camera_geometry)

add_executable(test_off_axis_orientation
    apps/test_off_axis_orientation.cpp
)
target_link_libraries(test_off_axis_orientation PRIVATE simlite_core)
add_test(NAME test_off_axis_orientation
         COMMAND test_off_axis_orientation ${CMAKE_CURRENT_SOURCE_DIR}/configs/mirror_1229_facets.csv)

add_executable(test_coordinate_frames
    apps/test_coordinate_frames.cpp
)
target_link_libraries(test_coordinate_frames PRIVATE simlite_core)
add_test(NAME test_coordinate_frames COMMAND test_coordinate_frames)

add_executable(test_light_collector_response
    apps/test_light_collector_response.cpp
)
target_link_libraries(test_light_collector_response PRIVATE simlite_core)
add_test(NAME test_light_collector_response COMMAND test_light_collector_response)

add_executable(test_detector_response_config
    apps/test_detector_response_config.cpp
)
target_link_libraries(test_detector_response_config PRIVATE simlite_core)
add_test(NAME test_detector_response_config COMMAND test_detector_response_config)

add_executable(test_telescope_optics
    apps/test_telescope_optics.cpp
)
target_link_libraries(test_telescope_optics PRIVATE simlite_core)
add_test(NAME test_telescope_optics COMMAND test_telescope_optics)

add_executable(test_photon_response_sampler
    apps/test_photon_response_sampler.cpp
)
target_link_libraries(test_photon_response_sampler PRIVATE simlite_core)
add_test(NAME test_photon_response_sampler COMMAND test_photon_response_sampler)

add_executable(test_trigger_response
    apps/test_trigger_response.cpp
)
target_link_libraries(test_trigger_response PRIVATE simlite_core)
add_test(NAME test_trigger_response COMMAND test_trigger_response)

add_executable(test_detector_pipeline
    apps/test_detector_pipeline.cpp
)
target_link_libraries(test_detector_pipeline PRIVATE simlite_core)
add_test(NAME test_detector_pipeline COMMAND test_detector_pipeline)

add_executable(test_array_timing_correction
    apps/test_array_timing_correction.cpp
)
target_link_libraries(test_array_timing_correction PRIVATE simlite_core)
add_test(NAME test_array_timing_correction COMMAND test_array_timing_correction)

add_executable(test_efficiency_curves
    apps/test_efficiency_curves.cpp
)
target_link_libraries(test_efficiency_curves PRIVATE simlite_core)
add_test(NAME test_efficiency_curves COMMAND test_efficiency_curves)

add_test(NAME compute_nsb_rate_spectral
         COMMAND compute_nsb_rate ${CMAKE_CURRENT_SOURCE_DIR}/configs/nsb/spectral_rate_check_with_obstruction.cfg)
set_tests_properties(compute_nsb_rate_spectral PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_atmosphere_transmission
    apps/test_atmosphere_transmission.cpp
)
target_link_libraries(test_atmosphere_transmission PRIVATE simlite_core)
add_test(NAME test_atmosphere_transmission
         COMMAND test_atmosphere_transmission
         ${CMAKE_CURRENT_SOURCE_DIR}/configs/atmosphere/atm_trans_4400_1_10_0_0_4400.dat)

add_executable(test_atmosphere_absorption_distribution
    apps/test_atmosphere_absorption_distribution.cpp
)
target_link_libraries(test_atmosphere_absorption_distribution PRIVATE simlite_core)
add_test(NAME test_atmosphere_absorption_distribution
         COMMAND test_atmosphere_absorption_distribution
         ${CMAKE_CURRENT_SOURCE_DIR}/configs/atmosphere/atm_trans_4400_1_10_0_0_4400.dat
         ${CMAKE_CURRENT_BINARY_DIR}/test_outputs/atmosphere_absorption_distribution.csv)

add_executable(scan_light_collector_angular_response
    apps/scan_light_collector_angular_response.cpp
)
target_link_libraries(scan_light_collector_angular_response PRIVATE simlite_core)

add_test(NAME run_optical_sim_generated
         COMMAND run_optical_sim ${CMAKE_CURRENT_SOURCE_DIR}/configs/test_generated.cfg)
set_tests_properties(run_optical_sim_generated PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME run_optical_sim_mirror_1229
         COMMAND run_optical_sim ${CMAKE_CURRENT_SOURCE_DIR}/configs/mirror_1229.cfg)
set_tests_properties(run_optical_sim_mirror_1229 PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME run_optical_sim_photon_csv
         COMMAND run_optical_sim ${CMAKE_CURRENT_SOURCE_DIR}/configs/mirror_1229_photon_csv_test.cfg)
set_tests_properties(run_optical_sim_photon_csv PROPERTIES
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(NAME run_optical_sim_photon_csv_filter
         COMMAND run_optical_sim ${CMAKE_CURRENT_SOURCE_DIR}/configs/mirror_1229_photon_csv_filter_test.cfg)
set_tests_properties(run_optical_sim_photon_csv_filter PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_test(
    NAME run_optical_sim_eventio_2d_obstruction
    COMMAND run_optical_sim
            ${CMAKE_CURRENT_SOURCE_DIR}/configs/test_eventio_2d_obstruction.cfg)
set_tests_properties(run_optical_sim_eventio_2d_obstruction PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    PASS_REGULAR_EXPRESSION "blocked_incoming=1")

add_test(NAME run_optical_sim_rejects_eventio
         COMMAND run_optical_sim ${CMAKE_CURRENT_SOURCE_DIR}/configs/official_tests/corsika_whiteboard.cfg)
set_tests_properties(run_optical_sim_rejects_eventio PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    WILL_FAIL TRUE)
