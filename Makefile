# Minimal user build wrapper.

BUILD_DIR ?= build
BUILD_TYPE ?= Release
HESSIO_ROOT ?= $(CURDIR)/external/hessioxxx/source
CMAKE ?= cmake
JOBS ?= $(shell command -v nproc >/dev/null 2>&1 && nproc || sysctl -n hw.ncpu 2>/dev/null || echo 4)
CMAKE_OSX_ARCHITECTURES ?= $(shell uname -m)
LACT_ENABLE_HDF5 ?= ON
LACT_ENABLE_ROOT ?= ON

.PHONY: all hessio hdf5 configure build clean distclean no-hessio

all: build

hessio:
	bash tools/build_hessio.sh

hdf5:
ifeq ($(LACT_ENABLE_HDF5),ON)
	bash tools/ensure_hdf5.sh "$(CMAKE_OSX_ARCHITECTURES)" "$(CMAKE)"
endif

configure: hessio hdf5
	$(CMAKE) -S . -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" -DLACT_ENABLE_HESSIO=ON -DLACT_ENABLE_HDF5="$(LACT_ENABLE_HDF5)" -DLACT_ENABLE_ROOT="$(LACT_ENABLE_ROOT)" -DHESSIO_ROOT="$(HESSIO_ROOT)" -DCMAKE_OSX_ARCHITECTURES="$(CMAKE_OSX_ARCHITECTURES)"

build: configure
	$(CMAKE) --build "$(BUILD_DIR)" -j$(JOBS)

no-hessio:
	$(CMAKE) -S . -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" -DLACT_ENABLE_HESSIO=OFF -DLACT_ENABLE_HDF5=OFF -DLACT_ENABLE_ROOT=OFF -DCMAKE_OSX_ARCHITECTURES="$(CMAKE_OSX_ARCHITECTURES)"
	$(CMAKE) --build "$(BUILD_DIR)" -j$(JOBS)

clean:
	$(CMAKE) --build "$(BUILD_DIR)" --target clean || true
	rm -rf "$(BUILD_DIR)"

distclean: clean
	rm -rf "$(HESSIO_ROOT)/out" "$(HESSIO_ROOT)/bin" "$(HESSIO_ROOT)/lib"
