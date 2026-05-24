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
    if ! "$CMAKE_BIN" -S "$TMP_DIR" -B "$TMP_DIR/build" -DCMAKE_OSX_ARCHITECTURES="$ARCH" >/tmp/lact_hdf5_check.log 2>&1; then
        return 1
    fi

    if [[ "$(uname -s)" == "Darwin" ]]; then
        local lib
        lib="$(grep -E '^HDF5_.*LIBRARY[^=]*:FILEPATH=' "$TMP_DIR/build/CMakeCache.txt" | head -1 | cut -d= -f2- || true)"
        if [[ -n "$lib" && -f "$lib" ]]; then
            if ! file "$lib" | grep -q "$ARCH"; then
                echo "Found HDF5 library, but its architecture does not include $ARCH:" >&2
                file "$lib" >&2
                return 2
            fi
        fi
    fi
}

install_hdf5() {
    case "$(uname -s)" in
        Darwin)
            local brew_bin
            if [[ "$ARCH" == "arm64" && -x /opt/homebrew/bin/brew ]]; then
                brew_bin=/opt/homebrew/bin/brew
            elif [[ "$ARCH" == "x86_64" && -x /usr/local/bin/brew ]]; then
                brew_bin=/usr/local/bin/brew
            elif command -v brew >/dev/null 2>&1; then
                brew_bin="$(command -v brew)"
            else
                brew_bin=""
            fi
            if [[ -n "$brew_bin" ]]; then
                "$brew_bin" install hdf5
            else
                echo "HDF5 not found and Homebrew is unavailable. Please install HDF5 C development libraries." >&2
                return 1
            fi
            ;;
        Linux)
            if command -v apt-get >/dev/null 2>&1; then
                sudo apt-get update
                sudo apt-get install -y libhdf5-dev
            elif command -v dnf >/dev/null 2>&1; then
                sudo dnf install -y hdf5-devel
            elif command -v yum >/dev/null 2>&1; then
                sudo yum install -y hdf5-devel
            else
                echo "HDF5 not found and no supported package manager was detected." >&2
                return 1
            fi
            ;;
        *)
            echo "HDF5 not found on this platform. Please install HDF5 C development libraries." >&2
            return 1
            ;;
    esac
}

if check_hdf5; then
    exit 0
fi

echo "HDF5 C library was not found for architecture $ARCH. Trying to install it..." >&2
if ! install_hdf5; then
    echo "Automatic HDF5 installation failed. Check download/install permissions, or build with LACT_ENABLE_HDF5=OFF." >&2
    exit 1
fi

if ! check_hdf5; then
    echo "HDF5 is still unavailable after the install attempt. Check architecture and permissions, or build with LACT_ENABLE_HDF5=OFF." >&2
    exit 1
fi
