
# PSC

[![Documentation Status](https://readthedocs.org/projects/psc/badge/?version=latest)](https://psc.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.com/psc-code/psc.svg?branch=master)](https://travis-ci.com/psc-code/psc)

## Introduction

This is the development version of PSC. It is usable for certain
problems, but overall lots of changes are happening and things may be
unstable / evolving.

PSC is a 3-dimensional fully electromagnetic particle-in-cell code for
kinetic plasma simulations. It supports Nvidia GPUs and patch-based
dynamic load balancing.

## Quickstart

### Building the code

PSC now uses cmake as its buildsystem. It is strongly recommended to build the code outside of the source in a different build directory. This is what I usually do:

```sh
$ cd src/psc # This is where the cloned git repository is (ie, the source)
$ mkdir build # This is my build directory
$ # create a cmake.sh so I don't have to remember the options I used
$ cat > cmake.sh <<EOF
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DUSE_CUDA=OFF \
    -DUSE_VPIC=OFF \
    ..
EOF
$ . cmake.sh # run cmake (hope for the best)
$ make

```

The above should build a lot of tests (you can use `ctest .` to run them) and two actual cases, `psc_bubble_yz` and `psc_flatfoil_yz` which you'll find in the `src/` subdirectory inside your build directory.

### Running a case

Make a directory for your run and change into it. This generation of PSC uses (almost) no command line options, rather everything gets hardcoded into a case (a.k.a. "input deck" in VPIC speak), e.g., `psc_flatfoil_yz.cxx`. (Eventually, the build system should be set up where a standalone case can be easily compiled, rather than expecting for it to be in the PSC source directory). So it should be as easy as 
```
mpirun -n 4 path/to/build/src/psc_flatfoil_yz
```

if you don't have to go through a batch system.

## Documentation

Well, this is really a big todo, so for now a lot of it is probably
all about emailing kai.germaschewski@unh.edu. But I'm working on
actually usable documentation and there's a start here: https://psc.readthedocs.io/en/latest/
