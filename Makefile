# Convenience wrapper for the minimal user build.
#
# Common usage:
#   make          # build hessioxxx and LACT_sim in ./build
#   make clean    # remove CMake build output only
#   make distclean# remove CMake build output and hessioxxx build products

BUILD_DIR ?= build
BUILD_TYPE ?= Release
HESSIO_ROOT ?= $(CURDIR)/external/hessioxxx/source
CMAKE ?= cmake
JOBS ?= $(shell command -v nproc >/dev/null 2>&1 && nproc || sysctl -n hw.ncpu 2>/dev/null || echo 4)
CMAKE_OSX_ARCHITECTURES ?= $(shell uname -m)
LACT_ENABLE_HDF5 ?= ON

.PHONY: all hessio hdf5 configure build clean distclean no-hessio

all: build

hessio:
	tools/build_hessio.sh

hdf5:
ifeq ($(LACT_ENABLE_HDF5),ON)
	tools/ensure_hdf5.sh "$(CMAKE_OSX_ARCHITECTURES)" "$(CMAKE)"
endif

configure: hessio hdf5
	$(CMAKE) -S . -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" -DLACT_ENABLE_HESSIO=ON -DLACT_ENABLE_HDF5="$(LACT_ENABLE_HDF5)" -DHESSIO_ROOT="$(HESSIO_ROOT)" -DCMAKE_OSX_ARCHITECTURES="$(CMAKE_OSX_ARCHITECTURES)"

build: configure
	$(CMAKE) --build "$(BUILD_DIR)" -j$(JOBS)

no-hessio:
	$(CMAKE) -S . -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" -DLACT_ENABLE_HESSIO=OFF -DCMAKE_OSX_ARCHITECTURES="$(CMAKE_OSX_ARCHITECTURES)"
	$(CMAKE) --build "$(BUILD_DIR)" -j$(JOBS)

clean:
	$(CMAKE) --build "$(BUILD_DIR)" --target clean || true
	rm -rf "$(BUILD_DIR)"

distclean: clean
	rm -rf "$(HESSIO_ROOT)/out" "$(HESSIO_ROOT)/bin" "$(HESSIO_ROOT)/lib"
