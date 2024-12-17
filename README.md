# OpenCafeMol

![Build Status](https://github.com/yutakasi634/OpenCafeMol/actions/workflows/main.yml/badge.svg)

## Requirement
- OpenMM
- C++17
- CMake 3.10

## Build
Since OpenCafeMol manages depending library via git submodule, clone this repo using git and do not download zip or release-tarball. It will cause compilation error.

```sh
$ git clone --recurse-submodules https://github.com/yutakasi634/OpenCafeMol.git
$ cd OpenCafeMol
$ mkdir build
$ cd build
$ cmake ..
$ make
```
After this, you will find an excutable binary in `bin` directory.  
If your OpenMM is not installed under `/usr` or `/usr/local`, you need to specify the install directory like `cmake .. -DOPENMM_ROOT=<your/openmm/installed/path>` in cmake step.

## Example
Example input files is in the sample directory. You can run this with the command below.
```sh
$ ./bin/open_cafemol sample/toml_input/sh3_AICG2+.toml # run simulation using toml interface
$ ./bin/open_cafemol sample/genesis_input/1SRL_cg.inp  # run simulation using genesis interface
```
The result files will be generated in output directory.
