# OpenAICG2plus

![Build Status](https://github.com/yutakasi634/OpenAICG2plus/actions/workflows/main.yml/badge.svg)

## Requirement
- OpenMM

## Build
Since OpenAICG2plus manages depending library via git submodule, clone this repo using git and do not download zip or release-tarball. It will cause compilation error.

```sh
$ git clone --recurse-submodules https://github.com/yutakasi634/OpenAICG2plus.git
$ cd OpenAICG2plus
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
$ mkdir output
$ ./bin/open_aicg2plus sample/toml_input/sh3_AICG2+.toml # run simulation using toml interface
$ ./bin/open_aicg2plus sample/genesis_input/1SRL_cg.inp  # run simulation using genesis interface
```
The result files will be generated in output directory.
