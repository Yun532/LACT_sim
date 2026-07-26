#!/usr/bin/env bash
set -euo pipefail

ARCH="${1:-$(uname -m)}"
CMAKE_BIN="${2:-cmake}"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

cat > "$TMP_DIR/CMakeLists.txt" <<'CMAKE'
cmake_minimum_required(VERSION 3.16)
project(check_hdf5 LANGUAGES C)
find_package(HDF5 COMPONENTS C REQUIRED)
CMAKE

check_hdf5() {
    rm -rf "$TMP_DIR/build"
    "$CMAKE_BIN" -S "$TMP_DIR" -B "$TMP_DIR/build" \
        -DCMAKE_OSX_ARCHITECTURES="$ARCH" >/tmp/lact_hdf5_check.log 2>&1
}

if check_hdf5; then
    exit 0
fi

echo "HDF5 C library not found; trying to install it..." >&2
case "$(uname -s)" in
    Darwin)
        command -v brew >/dev/null 2>&1 && brew install hdf5
        ;;
    Linux)
        if command -v apt-get >/dev/null 2>&1; then
            sudo apt-get update && sudo apt-get install -y libhdf5-dev
        elif command -v dnf >/dev/null 2>&1; then
            sudo dnf install -y hdf5-devel
        elif command -v yum >/dev/null 2>&1; then
            sudo yum install -y hdf5-devel
        else
            exit 1
        fi
        ;;
    *) exit 1 ;;
esac

check_hdf5 || {
    echo "HDF5 is still unavailable; install its C development library or build with LACT_ENABLE_HDF5=OFF." >&2
    exit 1
}
