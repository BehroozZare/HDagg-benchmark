
![APM](https://badgen.net/github/license/micromatch/micromatch)
![example workflow](https://github.com/sympiler/lbc/actions/workflows/cmakeUbuntu.yml/badge.svg)
![example workflow](https://github.com/sympiler/lbc/actions/workflows/cmakeMac.yml/badge.svg)


# HDagg
 [Hybrid Aggregation of Loop-carried
Dependence Iterations in Sparse Matrix
Computations (HDagg)](https://www.cs.toronto.edu/~mmehride/papers/HDag.pdf). is a DAG partitioning/scheduling algorithm used for making sparse matrix loops parallel.
 It can be used within code generators or libraries. It is integrated into Sympiler framework.
 This repository is the opensource reference implementation of the IPDPS 2022 paper.
 For more information see [Sympiler documents](https://www.sympiler.com/docs/lbc/).

## Files

* `src/`: source code
* `cmake/` and `CMakeLists.txt`: CMake files
* `input/`: input folder where matrices reside
* `output/`: output folder to store data
* `demo/`: example folder to show HDagg usage
* `scripts/`: Python and Bash scripts for generating and processing results

## Install

### Prerequisites
First following items should be installed:
* CMake
* C++ compiler (GCC, ICC, or CLang)
* METIS (optional) dependency for running the demo efficiently
  and is handled by the cmake. If you have installed the package using
  a packet manager (e.g., apt of homebrew), CMake should be able to detect it.
  Otherwise, it installs METIS from source internally.
* OpenMP (optional) for running some parts of the code in parallel. If you
  use GCC/ICC then OpenMP should be supported natively. If you use Apple CLang,
  you probably need to install OpenMP using `homebrew install libomp`. You can
* also install LLVM usng `brew install llvm` which support OpenMP natively.

#### MKL
The cmake currently supports mkl library inside oneapi. To use it, please provide `$MKLROOT` environmental variable using `export MKLROOT=<your-address>/make/latest`.  
* For example: `export MKLROOT=/opt/intel/oneapi/mkl/latest/` 

#### Relative Works:

You can switch the `HDAGG_WITH_SPMP` and `HDAGG_WITH_DAGP` options in `CMakeLists.txt` 
to add the [SpMP](https://github.com/IntelLabs/SpMP) or [DAGP](https://github.com/GT-TDAlab/dagP). To use these tools,
after installation, provide their addresses using `$SPMPROOT` and `$DAGPROOT`.
* For example: `export SPMPROOT=/home/behrooz/SpMP`

### Build
Then build HDagg, using the following:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```


You can always set `-DCMAKE_CXX_COMPILER=` and `-DCMAKE_C_COMPILER=` to use
a different compiler. For example:
`cmake -DCMAKE_CXX_COMPILER=/usr/local/Cellar/gcc\@9/9.3.0_2/bin/g++-9
-DCMAKE_C_COMPILER=/usr/local/Cellar/gcc\@9/9.3.0_2/bin/gcc-9 ..`

## Run

For each executable file the inputs is as follows: `./<executable_address> <matrix_address> <num_threads>` 
for example: `./build/demo/HDAGG_SpTRSV input/apache2.mtx 10`

