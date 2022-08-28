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

__global__ void CopyBoxForceGPUKernel(int molCount,
                            int molsInBox0,
                            int box,
                            uint* gpu_molLookup,
                            uint* gpu_molBoxCount,
                            int* gpu_startAtomIdx,
                            double *gpu_aForcex_old,
                            double *gpu_aForcey_old,
                            double *gpu_aForcez_old,
                            double *gpu_mForcex_old,
                            double *gpu_mForcey_old,
                            double *gpu_mForcez_old,
                            double *gpu_mTorquex_old,
                            double *gpu_mTorquey_old,
                            double *gpu_mTorquez_old,
                            double *gpu_aForcex_new,
                            double *gpu_aForcey_new,
                            double *gpu_aForcez_new,
                            double *gpu_mForcex_new,
                            double *gpu_mForcey_new,
                            double *gpu_mForcez_new,
                            double *gpu_mTorquex_new,
                            double *gpu_mTorquey_new,
                            double *gpu_mTorquez_new){

    int threadID = blockIdx.x * blockDim.x + threadIdx.x;
    // Optimal alu usage/memory latency will probably be 1 warp/molecule
    uint warpIdx = threadID / warp_size;
    uint laneIdx = threadID % warp_size;
    if (warpIdx >= molCount)
        return;

    int molIndex = gpu_molLookup[box*molsInBox0 + warpIdx];
    gpu_mForcex_new[molIndex] = gpu_mForcex_old[molIndex];
    gpu_mForcey_new[molIndex] = gpu_mForcey_old[molIndex];
    gpu_mForcez_new[molIndex] = gpu_mForcez_old[molIndex];
    gpu_mTorquex_new[molIndex] = gpu_mTorquex_old[molIndex];
    gpu_mTorquey_new[molIndex] = gpu_mTorquey_old[molIndex];
    gpu_mTorquez_new[molIndex] = gpu_mTorquez_old[molIndex];
    for (int particleIndex = gpu_startAtomIdx[molIndex] + laneIdx; particleIndex < gpu_startAtomIdx[molIndex + 1]; particleIndex += warp_size ){
      gpu_aForcex_new[particleIndex] = gpu_aForcex_old[particleIndex];
      gpu_aForcey_new[particleIndex] = gpu_aForcey_old[particleIndex];
      gpu_aForcez_new[particleIndex] = gpu_aForcez_old[particleIndex];
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


void CUDAMemoryUtils::CallCopyBoxForces(VariablesCUDA *vars,
                     XYZArray const &coords,
                     const MoleculeLookup& molLookup,
                     uint const buffer_index,
                     uint const box){

  BufferAccess<DeviceArray<double>, double, buffers> aFx_old(*(vars->gpu_aFx), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> aFy_old(*(vars->gpu_aFy), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> aFz_old(*(vars->gpu_aFz), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mFx_old(*(vars->gpu_mFx), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mFy_old(*(vars->gpu_mFy), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mFz_old(*(vars->gpu_mFz), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mTx_old(*(vars->gpu_mTx), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mTy_old(*(vars->gpu_mTy), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mTz_old(*(vars->gpu_mTz), buffer_index);

  BufferAccess<DeviceArray<double>, double, buffers> aFx_new(*(vars->gpu_aFx), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> aFy_new(*(vars->gpu_aFy), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> aFz_new(*(vars->gpu_aFz), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mFx_new(*(vars->gpu_mFx), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mFy_new(*(vars->gpu_mFy), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mFz_new(*(vars->gpu_mFz), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mTx_new(*(vars->gpu_mTx), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mTy_new(*(vars->gpu_mTy), buffer_index);
  BufferAccess<DeviceArray<double>, double, buffers> mTz_new(*(vars->gpu_mTz), buffer_index);




  int atomCount = coords.Count();
  int molCount = molLookup.NumInBox(box);

  // Run the kernel
  int threadsPerBlock = 256;
  int blocksPerGrid = (int)((molCount * warp_size) / threadsPerBlock) + 1;

  CopyBoxForceGPUKernel<<< blocksPerGrid, threadsPerBlock>>>(molCount, 
                            molLookup.NumInBox(0),
                            box,
                            molLookup.molLookupGPU->GetMolLookup(),
                            molLookup.molLookupGPU->GetNumMolsInBox(),
                            molLookup.molLookupGPU->GetStartAtomIdx(),
                            aFx_old->get(),
                            aFy_old->get(),
                            aFz_old->get(),
                            mFx_old->get(),
                            mFy_old->get(),
                            mFz_old->get(),
                            mTx_old->get(),
                            mTy_old->get(),
                            mTz_old->get(),
                            aFx_new->get(),
                            aFy_new->get(),
                            aFz_new->get(),
                            mFx_new->get(),
                            mFy_new->get(),
                            mFz_new->get(),
                            mTx_new->get(),
                            mTy_new->get(),
                            mTz_new->get());

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


