# Codex build notes

These notes record the local environment used while aligning LACT_sim ROOT
output with pylast. They are meant to avoid rediscovering the same setup details.

## Repositories

- LACT_sim: `/Users/yun/Downloads/LACT_sim`
- pylast: `/Users/yun/Downloads/LACT_sim/external/pylast`

The two repositories are independent. Do not commit `external/pylast` through
the LACT_sim repository.

## Local ROOT and pylast runtime

The usable ROOT-backed pylast test runtime found on this machine is:

```bash
/private/tmp/lact-root-test/bin/python
/private/tmp/lact-root-test/bin/root-config
```

Useful checks:

```bash
/private/tmp/lact-root-test/bin/python -c "import pylast; print(pylast.__file__)"
/private/tmp/lact-root-test/bin/root-config --cflags --libs
```

The pylast `.venv-py39` currently cannot import pylast without fixing its
runtime library path:

```text
Library not loaded: @rpath/libCore.so
```

## Local build limitation

The local LACT_sim and pylast CMake builds currently hit architecture-mismatch
link errors on this Mac: x86_64 targets try to link arm64 hessio/ROOT libraries.
This is a local build-environment issue, not a syntax error in the edited code.

Use syntax-only checks for ROOT-dependent files when this happens:

```bash
c++ -std=c++20 -fsyntax-only \
  -Iinclude \
  -I/private/tmp/lact-root-test/include \
  src/io/LactEventRootWriter.cpp

c++ -std=c++20 -fsyntax-only \
  -Iexternal/pylast/include \
  -Iexternal/pylast/include/external \
  -Iexternal/pylast/root/include \
  -I/private/tmp/lact-root-test/include \
  external/pylast/root/LactEventSource.cpp
```

Also run Python syntax checks for pylast wrappers:

```bash
python3 -m py_compile \
  external/pylast/src/pylast/io/LactEventSource.py \
  external/pylast/src/pylast/visualize/event_visualizer.py
```

## Server expectation

On the server, use a consistent architecture toolchain and a ROOT-enabled build.
The ROOT writer changes should be compiled there with the normal ROOT-enabled
LACT_sim build, then validated by reading the produced `lact_events.root` with
pylast `LactEventSource`.
