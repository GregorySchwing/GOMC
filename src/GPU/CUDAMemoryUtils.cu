#include "CUDAMemoryUtils.cuh"


// Selectively clear the force/torque of a single box
__global__ void ZeroBoxForceGPUKernel(int molCount,
                            int molsInBox0,
                            int box,
                            uint* gpu_molLookup,
                            uint* gpu_molBoxCount,
                            int* gpu_startAtomIdx,
                            double *gpu_aForcex,
                            double *gpu_aForcey,
                            double *gpu_aForcez,
                            double *gpu_mForcex,
                            double *gpu_mForcey,
                            double *gpu_mForcez,
                            double *gpu_mTorquex,
                            double *gpu_mTorquey,
                            double *gpu_mTorquez){

    int threadID = blockIdx.x * blockDim.x + threadIdx.x;
    // Optimal alu usage/memory latency will probably be 1 warp/molecule
    uint warpIdx = threadID / warp_size;
    uint laneIdx = threadID % warp_size;
    if (warpIdx >= molCount)
        return;

    int molIndex = gpu_molLookup[box*molsInBox0 + warpIdx];
    gpu_mForcex[molIndex] = 0.0;
    gpu_mForcey[molIndex] = 0.0;
    gpu_mForcez[molIndex] = 0.0;
    gpu_mTorquex[molIndex] = 0.0;
    gpu_mTorquey[molIndex] = 0.0;
    gpu_mTorquez[molIndex] = 0.0;
    for (int particleIndex = gpu_startAtomIdx[molIndex] + laneIdx; particleIndex < gpu_startAtomIdx[molIndex + 1]; particleIndex += warp_size ){
      gpu_aForcex[particleIndex] = 0.0;
      gpu_aForcey[particleIndex] = 0.0;
      gpu_aForcez[particleIndex] = 0.0;
    }
}


void CUDAMemoryUtils::CallZeroBoxForces(VariablesCUDA *vars,
                     XYZArray const &coords,
                     const MoleculeLookup& molLookup,
                     uint const buffer_index,
                     uint const box){

  BufferAccess<DeviceArray<double>, double, buffers> aFx(*(vars->gpu_aFx), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> aFy(*(vars->gpu_aFy), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> aFz(*(vars->gpu_aFz), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mFx(*(vars->gpu_mFx), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mFy(*(vars->gpu_mFy), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mFz(*(vars->gpu_mFz), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mTx(*(vars->gpu_mTx), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mTy(*(vars->gpu_mTy), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mTz(*(vars->gpu_mTz), buffer_index);




  int atomCount = coords.Count();
  int molCount = molLookup.NumInBox(box);

  // Run the kernel
  int threadsPerBlock = 256;
  int blocksPerGrid = (int)((molCount * warp_size) / threadsPerBlock) + 1;

  ZeroBoxForceGPUKernel<<< blocksPerGrid, threadsPerBlock>>>(molCount, 
                            molLookup.NumInBox(0),
                            box,
                            molLookup.molLookupGPU->GetMolLookup(),
                            molLookup.molLookupGPU->GetNumMolsInBox(),
                            molLookup.molLookupGPU->GetStartAtomIdx(),
                            aFx->get(),
                            aFy->get(),
                            aFz->get(),
                            mFx->get(),
                            mFy->get(),
                            mFz->get(),
                            mTx->get(),
                            mTy->get(),
                            mTz->get());

}
// Clear the force/torque of both boxes
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


