# Convenience wrapper for normal server builds.
#
# Common usage:
#   make          # build hessioxxx and LACT_sim in ./build
#   make test     # run the compact SII regression check
#   make clean    # remove CMake build output only
#   make distclean# remove CMake build output and hessioxxx build products

BUILD_DIR ?= build
BUILD_TYPE ?= Release
HESSIO_ROOT ?= $(CURDIR)/external/hessioxxx/source
CMAKE ?= cmake
JOBS ?= $(shell command -v nproc >/dev/null 2>&1 && nproc || sysctl -n hw.ncpu 2>/dev/null || echo 4)

.PHONY: all hessio configure build test clean distclean no-hessio no-root

all: build

hessio:
	tools/build_hessio.sh

configure: hessio
	$(CMAKE) -S . -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" -DHESSIO_ROOT="$(HESSIO_ROOT)"

build: configure
	$(CMAKE) --build "$(BUILD_DIR)" -j$(JOBS)

test:
	python3 -m pytest -q tests/test_sii_unified.py

no-hessio:
	$(CMAKE) -S . -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" -DLACT_ENABLE_HESSIO=OFF
	$(CMAKE) --build "$(BUILD_DIR)" -j$(JOBS)

no-root: hessio
	$(CMAKE) -S . -B "$(BUILD_DIR)" -DCMAKE_BUILD_TYPE="$(BUILD_TYPE)" -DHESSIO_ROOT="$(HESSIO_ROOT)" -DLACT_ENABLE_ROOT=OFF
	$(CMAKE) --build "$(BUILD_DIR)" -j$(JOBS)

clean:
	$(CMAKE) --build "$(BUILD_DIR)" --target clean || true
	rm -rf "$(BUILD_DIR)"

distclean: clean
	rm -rf "$(HESSIO_ROOT)/out" "$(HESSIO_ROOT)/bin" "$(HESSIO_ROOT)/lib"
