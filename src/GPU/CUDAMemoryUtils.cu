#include "CUDAMemoryUtils.cuh"

void CUDAMemoryUtils::ZeroForces(VariablesCUDA *vars,
                     XYZArray const &coords,
                     const MoleculeLookup& molLookup,
                     uint const buffer_index){

  int atomCount = coords.Count();
  int molCount = molLookup.molLookupCount;

  BufferAccess<DeviceArray<double>, double, buffers> aFx(*(vars->gpu_aFx), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> aFy(*(vars->gpu_aFy), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> aFz(*(vars->gpu_aFz), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mFx(*(vars->gpu_mFx), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mFy(*(vars->gpu_mFy), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mFz(*(vars->gpu_mFz), buffer_index);


  cudaMemset(aFx->get(), 0.0, atomCount * sizeof(double));
  cudaMemset(aFy->get(), 0.0, atomCount * sizeof(double));
  cudaMemset(aFz->get(), 0.0, atomCount * sizeof(double));
  cudaMemset(mFx->get(), 0.0, molCount * sizeof(double));
  cudaMemset(mFy->get(), 0.0, molCount * sizeof(double));
  cudaMemset(mFz->get(), 0.0, molCount * sizeof(double));
}
