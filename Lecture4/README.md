# Lecture 4

This Lecture will briefly introduce the progresses of GPU acceleration of Perturbo, including data structures, benchmarks, and compilation. The real simulation with GPU acceperation will be introduced in Hands-on 4. 

## Data structures and benchmarks of GPU Perturbo

Please refer to the slides of Lecture 4.

## Compilation of GPU version of QE and PERTURBO
I will briefly introduce how to compile the GPU version of QE and PERTURBO on Perlmutter.

### Environment setup on Perlmutter
The environment set is very important for the compilation of GPU version of QE and PERTURBO. Please follow the steps below to set up the environment on Perlmutter. 

```bash
module load cudatoolkit/12.0 PrgEnv-nvidia nvidia/23.1
module load cray-hdf5-parallel/1.12.2.3 cray-libsci/23.02.1.1
module load conda python/3.11
module load cray-fftw/3.3.10.8
```

These instructions have been verified to work with NVHPC Toolkit 23.1. Specifically, we have identified an nvfortran compiler bug in version 23.9-24.5 that prevents us from using that version of the Toolkit.


### QE 
Before compiling PERTURBO, we need to compile a GPU version of Quantum ESPRESSO (QE) first. Instead of using docker for other sessions, you need download the source code from the github repository

Download the source code of QE from github and checkout to version 7.3.1,
```bash
cd ~/work/Perturbo-workshop-2025/Hands-on4
git clone https://github.com/QEF/q-e.git
cd q-e
git checkout qe-7.3.1
```

Move on the configure step,
```bash
FC=nvfortran F90=nvfortran MPIF90=mpif90 CC=cc \
  FFLAGS="-fast -Mlarge_arrays -mcmodel=medium" CFLAGS="-fast -mcmodel=medium" \
  LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/lib64" \
  BLAS_LIBS="-lblas" LAPACK_LIBS="-llapack" \
  ./configure --enable-parallel=yes --enable-openmp=yes \
  --with-cuda=/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/cuda \
  --with-cuda-cc=80 --with-cuda-runtime=12.0 \
  --with-hdf5=/opt/cray/pe/hdf5-parallel/1.12.2.3/nvidia/20.7
```

Please check that the hdf5 library has been correctly linked. Sometimes, the configure script does not work properly.

```bash
HDF5_LIBS = -L/opt/cray/pe/hdf5-parallel/1.12.2.3/nvidia/20.7/lib -lhdf5 -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl
```


Compile the code using 20 threads,
```bash
make -j 20 pw ph pp w90
```

You can check whether the executables are generated in the bin directory,
```bash
ls bin
```

### PERTURBO

```bash
cd ~/work/Perturbo-workshop-2025/Hands-on4/q-e 
git clone git@github.com:perturbo-code/perturbo.git
cd perturbo
cp config/make_nersc_nvhpc.sys ./make.sys
```
Besides changing the library path, two lines like below need to be uncommented, which should be the default setting.

```bash
FFLAGS += -acc=gpu -Minfo=accel -gpu=cc80,nomanaged,deepcopy,zeroinit
LDFLAGS += -acc=gpu -gpu=cc80,nomanaged,deepcopy,zeroinit
```
Finally, a simple command by executing
```bash
make
```

For more details about the flags in `make.sys`, please follow the presenter.

You can check whether the executables are generated in the bin directory,
```bash
ls bin
```
There should be `perturbo.x` and `qe2pert.x`

### Compilation flags for GPU codes
Except for the general flags in `make.sys`, there are some important flags specific for GPU codes of Perturbo. Please check them in `pert-src/pert_config.h`. Here are some introductions:

- `SCAT_FWD`: Defining `SCAT_FWD` causes PERTURBO to generate a data structure that supports forward traversal, i.e. "source element ordered" calculation. 
- `SCAT_REV`: Defining `SCAT_REV` causes PERTURBO to generate a data structure that supports reverse traversal, i.e. "target element ordered" calculation.
- `STORE_G2DT`: Defining `STORE_G2DT` causes the `scatter_channels` array to store the two `-g2*dt1` and `-g2*dt2` terms, instead of just holding g2.  This increases memory requirements of PERTURBO, but can result in faster computations.
- `CDYN_USE_ARRAY_REDUCE`: Defining `CDYN_USE_ARRAY_REDUCE` causes the OpenMP code to use an array-reduction instead of `omp atomic` to perform updates to the epcol array.  This means that each OpenMP thread has its own copy of epcol, which uses more memory, and combining them at the end will take time as well.  However, array-reduction is usually much faster than `omp atomic`.
- `CDYN_SORT_SCAT_TGTS`: Defining `CDYN_SORT_SCAT_TGTS` causes the target-oriented code to order the scatter-target array elements by decreasing sub-array lengths.  For the GPU this seems to greatly improve performance, since all threads in a warp will finish working around the same time as each other.

Note that the default values in _pert-src/pert\_config.h_ are suggested. Please try other options with care and make sure they align with your intended use. 

