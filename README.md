# Sparsity

Sparsity is a collection of benchmarks that can be used to evaluate the
performance of sparse linear algebra implementations. The code in this
repository is taken from open source benchmark suites, real world chemical
simulation software, and simulated workloads from other problem domains.

# Building Sparsity

The benchmarks in this package must be built and run separately using their own
build systems.

```
git clone https://github.com/Baltoli/sparsity.git
cd sparsity
```

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
SPMV_INSTALL_PREFIX=prefix make platform={native,mkl,gpu,opencl,clgpu} install
```

The installed library will be named `libplatform-spmv.so`.

## NPB

We use the conjugate gradient benchmark from the NAS benchmark suite. To build
and run it:

```
cd NPB3.3.1
```
