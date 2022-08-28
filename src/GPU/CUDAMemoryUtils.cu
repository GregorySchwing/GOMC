#include "CUDAMemoryUtils.cuh"

void CUDAMemoryUtils::ZeroForces(VariablesCUDA *vars,
                     XYZArray const &coords,
                     const MoleculeLookup& molLookup,
                     uint const buffer_index){

  int atomCount = coords.Count();
  int molCount = molLookup.molLookupCount;
  BufferAccess<DeviceArray<double>, double, buffers> gpu_LJEn(*(vars->gpu_LJEn), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> gpu_REn(*(vars->gpu_REn), buffer_index);

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
  cudaMemset(gpu_LJEn->get(), 0.0, 1 * sizeof(double));
  cudaMemset(gpu_REn->get(), 0.0, 1 * sizeof(double));

  BufferAccess<DeviceArray<double>, double, buffers> mTx(*(vars->gpu_mTx), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mTy(*(vars->gpu_mTy), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mTz(*(vars->gpu_mTz), buffer_index);

  cudaMemset(mTx->get(), 0.0, molCount * sizeof(double));
  cudaMemset(mTy->get(), 0.0, molCount * sizeof(double));
  cudaMemset(mTz->get(), 0.0, molCount * sizeof(double));
}
