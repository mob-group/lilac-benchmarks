###############################################################################
#
# CMake file for generating CUDA routines
#
###############################################################################

cmake_minimum_required(VERSION 2.8)
project(CUDAINTERFACE)

# Default to Release build
set(CMAKE_BUILD_TYPE Release)

# svn root directory should be one directory above this directory
get_filename_component(SVN_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/.. ABSOLUTE)

# CUDA compilation
find_package(CUDA)
# Set nvcc flags
# These use optimisations and build for CUDA 3.5 virtual and real architectures
# See http://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-compilation 
# for more details.
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-O3;-arch=compute_35;-code=sm_35)

# Make sure use_DPDP flag is set for AMBER
add_definitions(-Duse_DPDP)

# Set source and include directories
file(GLOB CUDA_SOURCES "${CUDAINTERFACE_SOURCE_DIR}/bfgsts/*.cu" "${CUDAINTERFACE_SOURCE_DIR}/cost_function/*.cu" "${CUDAINTERFACE_SOURCE_DIR}/lbfgs/*.cu" "${CUDAINTERFACE_SOURCE_DIR}/potential/*.cu" "${CUDAINTERFACE_SOURCE_DIR}/utils/*.cu")
set(INCLUDE_DIRS "${CUDAINTERFACE_SOURCE_DIR}/include" "${CUDAINTERFACE_SOURCE_DIR}/lbfgs" "${SVN_ROOT}/AMBER12/cuda")

# Create library
cuda_include_directories(${INCLUDE_DIRS})
cuda_add_library(CUDAinterface ${CUDA_SOURCES})

# Link against cuBLAS
cuda_add_cublas_to_target(CUDAinterface)
