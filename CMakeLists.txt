cmake_minimum_required(VERSION 3.16)
project(simlite LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(include)
enable_testing()

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

add_library(simlite_core
    src/io/SyntheticPhotonSource.cpp
    src/io/PhotonCsvSource.cpp
    src/io/MirrorFacetCsvReader.cpp
    src/app/OpticalSimCommon.cpp
    src/geometry/FacetFactory.cpp
    src/optics/RayTracer.cpp
    src/optics/OpticalTracer.cpp
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

add_executable(run_optical_sim
    apps/run_optical_sim.cpp
)
target_link_libraries(run_optical_sim PRIVATE simlite_core)

if(LACT_HAS_HESSIO)
    add_executable(run_corsika_trace
        apps/run_corsika_trace.cpp
    )
    target_link_libraries(run_corsika_trace PRIVATE simlite_core)
    if(LACT_HAS_HDF5)
        target_compile_definitions(run_corsika_trace PRIVATE LACT_HAS_HDF5=1)
        target_include_directories(run_corsika_trace PRIVATE ${HDF5_INCLUDE_DIRS})
        target_link_libraries(run_corsika_trace PRIVATE ${HDF5_LIBRARIES})
    endif()
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
