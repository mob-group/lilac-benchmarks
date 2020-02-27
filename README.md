Benchmarks for evaluating the performance of sparse linear algebra
implementations. The code in this repository is taken from open source benchmark
suites, real world chemical simulation software, and simulated workloads from
other problem domains.

Results from these benchmarks appear in the following paper:

* [**Automatically Harnessing Sparse Acceleration**][cc] (CC 2020)

# Building

The benchmarks in this package must be built and run separately using their own
build systems.

```
git clone https://github.com/Baltoli/sparsity.git
cd sparsity
```

`LD_LIBRARY_PATH` needs to be able to find the installed `libspmv` libraries
(see below), as well as any required by the different platforms, both at compile
and run time.

## `libspmv`

To build the shared library linear algebra implementations:

```
cd libspmv
make {native,mkl,gpu,opencl,clgpu}.so
```

Each implementation requires different libraries to be installed. These
dependencies are:

* `native`: The sequential C implementation has no dependencies and should build
  on any system with a C compiler to act as a performance baseline.
* `mkl`: Requires an Intel CPU with MKL installed. The environment variable
  `$MKLROOT` should point to the installation of MKL, and can be defined by
  running the shell scripts supplied with MKL itself.
* `gpu`: Requires a CUDA-compatible Nvidia GPU and CUDA 9.1 installed, with
  installation root `$CUDA_ROOT`.
* `opencl`: OpenCL version that uses the first available OpenCL device on
  your platform, with OpenCL installed in `$CL_ROOT`. On our test platform, the
  first available device is the integrated GPU.
* `clgpu`: As for OpenCL, but using the second device - this corresponds to the
  external GPU on our test platform.

To install a built shared library:

```
SPMV_ROOT=prefix make platform={native,mkl,gpu,opencl,clgpu} install
```

The installed library will be named `libplatform-spmv.so`.

## NPB

We use the conjugate gradient benchmark from the NAS benchmark suite. To build
and run it:

```
cd NPB3.3.1
mkdir bin
make SPMV_VERSION=platform F77=gfortran CLASS=C CG
./bin/cg.C.x
```

The relevant benchmark results are the MOp/s and seconds taken to run. `F77` can
be set to whatever fortran compiler is most appropriate.

## SparseBench

This application implements sparse matrix benchmarks for various different
conjugate gradient variants. We use the general, unsymmetric version.

The first step is to generate suitably large matrices:
```
./big_gen.py --filename crsmat170u --size 170 write
```

[cc]: https://dl.acm.org/doi/10.1145/3377555.3377893
