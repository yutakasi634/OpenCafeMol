# OpenAICG2plus

## Requirement
- OpenMM

## Build
Since OpenAICG2plus manages depending library via git submodule, clone this repo using git and do not download zip or release-tarball. It will cause compilation error.

```sh
$ git clone https://github.com/yutakasi634/OpenAICG2plus.git
$ cd OpenAICG2plus
$ mkdir build
$ cd build
$ cmake ..
$ make
```
After this, you will find an excutable binary in `bin` directory.  
If your OpenMM is not install under `/usr` or `/usr/local`, you need to specify the install directory like `cmake .. -DOPRNMM_ROOT` in cmake step.
