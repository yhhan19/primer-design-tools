# Installation and build

[Quick start](../README.md)

## Build

Use Linux, macOS, or WSL with Git, GNU Make, a C++17 compiler, and POSIX threads.
Python 3 is needed for regression tests and the optional comparison plot, not the
main executable. Run from the repository root:

```bash
git submodule update --init --recursive
make -C primer3/src
make -j2
./bin/dpro --help
```

The root Makefile links Primer3 object files and does **not** build them
automatically. Build Primer3 first. The executable is `bin/dpro`.

Thermodynamic parameters come from `PRIMER3_PATH/src/primer3_config`. The default
build uses `./primer3`, so keep that directory available and run from the
repository root. To use another checkout, build it first, then use
`make PRIMER3=/absolute/path/to/primer3`; rebuild application objects when changing
this compiled-in path.
