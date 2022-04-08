# OpenAICG2plus

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
If your OpenMM is not installed under `/usr` or `/usr/local`, you need to specify the install directory like `cmake .. -DOPRNMM_ROOT=<your/openmm/installed/path>` in cmake step.

## Example
Example input files is in the input directory. You can run this with the command below.
```sh
$ ./bin/open_aicg2plus input/sh3_AICG2+.toml
```
The result files will be generated in output directory.
