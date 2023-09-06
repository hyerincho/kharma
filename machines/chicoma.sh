# LANL Machines: HPC and IC

# Chicoma
if [[ "$HOST" == "ch-fe"* ]]; then
  HOST_ARCH="ZEN2"

  # Cray environments get confused easy
  # Make things as simple as possible
  # TODO version with Cray wrappers?
  module purge
  export CRAY_CPU_TARGET="x86-64"
  if [[ "$ARGS" == *"cuda"* ]]; then
    DEVICE_ARCH="AMPERE80"
    # System HDF5 can't use compression
    EXTRA_FLAGS="-DPARTHENON_DISABLE_HDF5_COMPRESSION=ON $EXTRA_FLAGS"
    # Runtime
    MPI_NUM_PROCS=4
    if [[ "$ARGS" == *"gnu"* ]]; then
      module load PrgEnv-gnu cpe-cuda cuda
    elif [[ "$ARGS" == *"intel"* ]]; then
      module load PrgEnv-intel
    elif [[ "$ARGS" == *"nvc++"* ]]; then
      module load PrgEnv-nvhpc cray-hdf5-parallel
      EXTRA_FLAGS="-DCMAKE_CUDA_COMPILER=$HOME/bin/nvc++-wrapper -DCMAKE_CUDA_COMPILER_ID=NVHPC -DCMAKE_CUDA_COMPILER_VERSION=11.6 $EXTRA_FLAGS"
    else
      module load PrgEnv-nvhpc cray-hdf5-parallel
    fi
  else
    module load PrgEnv-aocc
  fi
  module load cmake

  # Runtime
  MPI_NUM_PROCS=4
  MPI_EXE=srun
  MPI_EXTRA_ARGS="--cpu-bind=mask_cpu:0x0*16,0x1*16,0x2*16,0x3*16 ~/bin/select-gpu"
  unset OMP_NUM_THREADS
  unset OMP_PROC_BIND
  unset OMP_PLACES
fi
