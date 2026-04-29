# Project-Local hessioxxx Dependency

This directory vendors the hessioxxx/EventIO library used by LACT_sim's
CORSIKA/EventIO photon-bunch converter.

## Source

Official upstream location:

```text
https://www.mpi-hd.mpg.de/hfm/bernlohr/sim_telarray/hessioxxx.tar.gz
https://www.mpi-hd.mpg.de/hfm/bernlohr/sim_telarray/sha256.sum
```

Downloaded archive:

```text
external/hessioxxx/hessioxxx.tar.gz
```

Official SHA-256:

```text
3c6e0221112a7d0bbc35ef08213c6cb28c3268fab462ca3703422dbc8acd815c  hessioxxx.tar.gz
```

Local verification:

```bash
shasum -a 256 external/hessioxxx/hessioxxx.tar.gz
```

The current archive's `ChangeLog` starts at:

```text
2025-02-20  Konrad Bernloehr
```

## Build

From the LACT_sim repository root:

```bash
tools/build_hessio.sh
```

The script creates the `out/`, `bin/`, and `lib/` directories before running
the upstream Makefile. If running hessioxxx manually, the equivalent is:

```bash
mkdir -p external/hessioxxx/source/out external/hessioxxx/source/bin external/hessioxxx/source/lib
make -C external/hessioxxx/source
```

The build produces:

```text
external/hessioxxx/source/bin/read_iact
external/hessioxxx/source/lib/libhessio.so    # Linux
external/hessioxxx/source/lib/libhessio.dylib # macOS
```

On macOS, runtime commands that link `libhessio.dylib` need:

```bash
export DYLD_LIBRARY_PATH=external/hessioxxx/source/lib
```

On Linux, use:

```bash
export LD_LIBRARY_PATH=external/hessioxxx/source/lib
```

## LACT_sim EventIO Converter

Build LACT_sim's converter against this vendored hessioxxx:

```bash
./tools/build_eventio_converter.sh
```

Convert CORSIKA/EventIO photon bunches:

```bash
DYLD_LIBRARY_PATH=external/hessioxxx/source/lib \
  ./build_eventio/eventio_to_photon_csv \
  input.zst \
  output.csv \
  --event-id-mode event_array100
```

For the validation file
`/Users/yun/Downloads/photon_E500_th0_run000001.zst`, the project-local build
wrote 1245159 photon-bunch rows and matched the provided ROOT conversion over
the full table.

Comparison log:

```text
run_logs/root_bunch_adapter/eventio_direct_external_vs_root_full_comparison.txt
```
