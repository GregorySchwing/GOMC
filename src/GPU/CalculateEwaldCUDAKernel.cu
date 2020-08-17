/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifdef GOMC_CUDA

#include <cuda.h>
#include <cuda_runtime.h>
#include "BoxDimensions.h"
#include "CalculateEwaldCUDAKernel.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "CalculateMinImageCUDAKernel.cuh"
#include "CUDAMemoryManager.cuh"
#include "cub/cub.cuh"
#include <vector>

using namespace cub;

#define FULL_MASK 0xffffffff

void CallBoxReciprocalSetupGPU(VariablesCUDA *vars,
                               XYZArray const &coords,
                               double const *kx,
                               double const *ky,
                               double const *kz,
                               std::vector<double> particleCharge,
                               uint imageSize,
                               double *sumRnew,
                               double *sumInew,
                               double *prefact,
                               double *hsqr,
                               double &energyRecip,
                               uint box)
{
  double *gpu_particleCharge;
  double * gpu_energyRecip;
  double * gpu_final_energyRecip;
  int blocksPerGrid, threadsPerBlock;
  int atomNumber = coords.Count();

  CUMALLOC((void**) &gpu_particleCharge,
           particleCharge.size() * sizeof(double));
  CUMALLOC((void**) &gpu_energyRecip, imageSize * sizeof(double));
  CUMALLOC((void**) &gpu_final_energyRecip, sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(double),
             cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_kx[box], kx, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_ky[box], ky, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_kz[box], kz, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_prefact[box], prefact, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);
  cudaMemcpy(vars->gpu_hsqr[box], hsqr, imageSize * sizeof(double),
             cudaMemcpyHostToDevice);
  checkLastErrorCUDA(__FILE__, __LINE__);

  threadsPerBlock = 256;
  blocksPerGrid = (int)(imageSize / threadsPerBlock) + 1;
  BoxReciprocalSetupGPU <<< blocksPerGrid, threadsPerBlock>>>(
    vars->gpu_x,
    vars->gpu_y,
    vars->gpu_z,
    vars->gpu_kx[box],
    vars->gpu_ky[box],
    vars->gpu_kz[box],
    atomNumber,
    gpu_particleCharge,
    vars->gpu_sumRnew[box],
    vars->gpu_sumInew[box],
    imageSize);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  BoxReciprocalGPU <<< blocksPerGrid, threadsPerBlock>>>(
    vars->gpu_prefact[box],
    vars->gpu_sumRnew[box],
    vars->gpu_sumInew[box],
    gpu_energyRecip,
    imageSize);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(sumRnew, vars->gpu_sumRnew[box],
             imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, vars->gpu_sumInew[box],
             imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);

  // ReduceSum
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecip,
                    gpu_final_energyRecip, imageSize);
  CUMALLOC(&d_temp_storage, temp_storage_bytes);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecip,
                    gpu_final_energyRecip, imageSize);
  cudaMemcpy(&energyRecip, gpu_final_energyRecip,
             sizeof(double), cudaMemcpyDeviceToHost);

  CUFREE(gpu_particleCharge);
  CUFREE(gpu_energyRecip);
  CUFREE(gpu_final_energyRecip);
  CUFREE(d_temp_storage);
}

void CallMolReciprocalGPU(VariablesCUDA *vars,
                          XYZArray const &currentCoords,
                          XYZArray const &newCoords,
                          std::vector<double> particleCharge,
                          uint imageSize,
                          double *sumRnew,
                          double *sumInew,
                          double &energyRecipNew,
                          uint box)
{
  // Calculate atom number
  int atomNumber = currentCoords.Count();
  int newCoordsNumber = newCoords.Count();
  double *gpu_particleCharge;
  int blocksPerGrid, threadsPerBlock;
  double *gpu_energyRecipNew, *gpu_final_energyRecipNew;

  CUMALLOC((void**) &gpu_particleCharge,
           particleCharge.size() * sizeof(double));
  CUMALLOC((void**) &gpu_energyRecipNew, imageSize * sizeof(double));
  CUMALLOC((void**) &gpu_final_energyRecipNew, sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(double),
             cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_x, currentCoords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, currentCoords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, currentCoords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_nx, newCoords.x, newCoordsNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_ny, newCoords.y, newCoordsNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_nz, newCoords.z, newCoordsNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  threadsPerBlock = 256;
  blocksPerGrid = (int)(imageSize / threadsPerBlock) + 1;
  MolReciprocalGPU <<< blocksPerGrid,
                   threadsPerBlock>>>(vars->gpu_x, vars->gpu_y, vars->gpu_z,
                                      vars->gpu_nx, vars->gpu_ny, vars->gpu_nz,
                                      vars->gpu_kxRef[box], vars->gpu_kyRef[box],
                                      vars->gpu_kzRef[box],
                                      atomNumber,
                                      gpu_particleCharge,
                                      vars->gpu_sumRnew[box],
                                      vars->gpu_sumInew[box],
                                      vars->gpu_sumRref[box],
                                      vars->gpu_sumIref[box],
                                      vars->gpu_prefactRef[box],
                                      gpu_energyRecipNew,
                                      imageSize);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(sumRnew, vars->gpu_sumRnew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, vars->gpu_sumInew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);

  // ReduceSum
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecipNew,
                    gpu_final_energyRecipNew, imageSize);
  CUMALLOC(&d_temp_storage, temp_storage_bytes);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecipNew,
                    gpu_final_energyRecipNew, imageSize);
  cudaMemcpy(&energyRecipNew, gpu_final_energyRecipNew,
             sizeof(double), cudaMemcpyDeviceToHost);


  CUFREE(gpu_particleCharge);
  CUFREE(gpu_energyRecipNew);
  CUFREE(gpu_final_energyRecipNew);
  CUFREE(d_temp_storage);
}

void CallSwapReciprocalGPU(VariablesCUDA *vars,
                           XYZArray const &coords,
                           std::vector<double> particleCharge,
                           uint imageSize,
                           double *sumRnew,
                           double *sumInew,
                           int const insert,
                           double &energyRecipNew,
                           uint box)
{
  // Calculate atom number
  int atomNumber = coords.Count();
  // given coordinates
  double *gpu_particleCharge;
  int blocksPerGrid, threadsPerBlock;
  double *gpu_energyRecipNew, *gpu_final_energyRecipNew;

  CUMALLOC((void**) &gpu_particleCharge,
           particleCharge.size() * sizeof(double));
  CUMALLOC((void**) &gpu_energyRecipNew, imageSize * sizeof(double));
  CUMALLOC((void**) &gpu_final_energyRecipNew, sizeof(double));

  cudaMemcpy(gpu_particleCharge, &particleCharge[0],
             particleCharge.size() * sizeof(double),
             cudaMemcpyHostToDevice);

  cudaMemcpy(vars->gpu_x, coords.x, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, coords.y, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, coords.z, atomNumber * sizeof(double),
             cudaMemcpyHostToDevice);

  threadsPerBlock = 256;
  blocksPerGrid = (int)(imageSize / threadsPerBlock) + 1;
  SwapReciprocalGPU <<< blocksPerGrid,
                    threadsPerBlock>>>(vars->gpu_x, vars->gpu_y, vars->gpu_z,
                                       vars->gpu_kxRef[box], vars->gpu_kyRef[box],
                                       vars->gpu_kzRef[box],
                                       atomNumber,
                                       gpu_particleCharge,
                                       vars->gpu_sumRnew[box],
                                       vars->gpu_sumInew[box],
                                       vars->gpu_sumRref[box],
                                       vars->gpu_sumIref[box],
                                       vars->gpu_prefactRef[box],
                                       insert,
                                       gpu_energyRecipNew,
                                       imageSize);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
//#ifndef NDEBUG
  // In the future maybe we could remove this for Nondebug?
  cudaMemcpy(sumRnew, vars->gpu_sumRnew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(sumInew, vars->gpu_sumInew[box], imageSize * sizeof(double),
             cudaMemcpyDeviceToHost);
//#endif

  // ReduceSum
  void *d_temp_storage = NULL;
  size_t temp_storage_bytes = 0;
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecipNew,
                    gpu_final_energyRecipNew, imageSize);
  CUMALLOC(&d_temp_storage, temp_storage_bytes);
  DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, gpu_energyRecipNew,
                    gpu_final_energyRecipNew, imageSize);
  cudaMemcpy(&energyRecipNew, gpu_final_energyRecipNew,
             sizeof(double), cudaMemcpyDeviceToHost);

  CUFREE(gpu_particleCharge);
  CUFREE(gpu_energyRecipNew);
  CUFREE(gpu_final_energyRecipNew);
  CUFREE(d_temp_storage);
}

void CallBoxForceReciprocalGPU(
  VariablesCUDA *vars,
  XYZArray &atomForceRec,
  XYZArray &molForceRec,
  const std::vector<double> &particleCharge,
  const std::vector<int> &particleMol,
  const std::vector<int> &particleKind,
  const std::vector<bool> &particleHasNoCharge,
  const std::vector<int> &startMol,
  const std::vector<int> &lengthMol,
  double alpha,
  double alphaSq,
  double qqFact,
  double constValue,
  uint imageSize,
  XYZArray const &molCoords,
  int boxStart,
  int boxEnd,
  BoxDimensions const &boxAxes, 
  int box
)
{
  int numberOfAtomsInsideBox = boxEnd - boxStart;
  int atomCount = atomForceRec.Count();
  int molCount = molForceRec.Count();
  double *gpu_particleCharge;
  int *gpu_particleMol, *gpu_particleKind;
  bool *gpu_particleHasNoCharge;
  bool *arr_particleHasNoCharge = new bool[particleHasNoCharge.size()];
  int *gpu_startMol, *gpu_lengthMol;

  // particleHasNoCharge is stored in vector<bool>, so in order to copy it to GPU
  // it needs to be stored in bool[]. because:
  // std::vector<bool> : Does not necessarily store its elements as a contiguous array
  for(int i=0; i<particleHasNoCharge.size(); i++) {
    arr_particleHasNoCharge[i] = particleHasNoCharge[i];
  }

  // calculate block and grid sizes
  int threadsPerBlock = 256;
  int blocksPerGrid = numberOfAtomsInsideBox;

  CUMALLOC((void **) &gpu_particleCharge, particleCharge.size() * sizeof(double));
  CUMALLOC((void **) &gpu_particleHasNoCharge, particleHasNoCharge.size() * sizeof(bool));
  CUMALLOC((void **) &gpu_startMol, startMol.size() * sizeof(int));
  CUMALLOC((void **) &gpu_lengthMol, lengthMol.size() * sizeof(int));
  CUMALLOC((void **) &gpu_particleMol, particleMol.size() * sizeof(int));
  CUMALLOC((void **) &gpu_particleKind, particleKind.size() * sizeof(int));

  cudaMemcpy(vars->gpu_aForceRecx, atomForceRec.x, sizeof(double) * atomCount, cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_aForceRecy, atomForceRec.y, sizeof(double) * atomCount, cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_aForceRecz, atomForceRec.z, sizeof(double) * atomCount, cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecx, molForceRec.x, sizeof(double) * molCount, cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecy, molForceRec.y, sizeof(double) * molCount, cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecz, molForceRec.z, sizeof(double) * molCount, cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleCharge, &particleCharge[0], sizeof(double) * particleCharge.size(), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleMol, &particleMol[0], sizeof(int) * particleMol.size(), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleKind, &particleKind[0], sizeof(int) * particleKind.size(), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleHasNoCharge, arr_particleHasNoCharge, sizeof(bool) * particleHasNoCharge.size(), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, molCoords.x, sizeof(double) * atomCount, cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, molCoords.y, sizeof(double) * atomCount, cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, molCoords.z, sizeof(double) * atomCount, cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_startMol, &startMol[0], sizeof(int) * startMol.size(), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_lengthMol, &lengthMol[0], sizeof(int) * lengthMol.size(), cudaMemcpyHostToDevice);

  checkLastErrorCUDA(__FILE__, __LINE__);
  BoxForceReciprocalGPU<<<blocksPerGrid, threadsPerBlock>>>(
    vars->gpu_aForceRecx,
    vars->gpu_aForceRecy,
    vars->gpu_aForceRecz,
    vars->gpu_mForceRecx,
    vars->gpu_mForceRecy,
    vars->gpu_mForceRecz,
    gpu_particleCharge,
    gpu_particleMol,
    gpu_particleKind,
    gpu_particleHasNoCharge,
    gpu_startMol,
    gpu_lengthMol,
    alpha,
    alphaSq,
    qqFact,
    constValue,
    imageSize,
    vars->gpu_kx[box],
    vars->gpu_ky[box],
    vars->gpu_kz[box],
    vars->gpu_x,
    vars->gpu_y,
    vars->gpu_z,
    vars->gpu_prefact[box],
    vars->gpu_sumRnew[box],
    vars->gpu_sumInew[box],
    vars->gpu_isFraction,
    vars->gpu_molIndex,
    vars->gpu_kindIndex,
    vars->gpu_lambdaCoulomb,
    boxAxes.GetAxis(box).x,
    boxAxes.GetAxis(box).y,
    boxAxes.GetAxis(box).z,
    box
  );
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(atomForceRec.x, vars->gpu_aForceRecx, sizeof(double) * atomCount, cudaMemcpyDeviceToHost);
  cudaMemcpy(atomForceRec.y, vars->gpu_aForceRecy, sizeof(double) * atomCount, cudaMemcpyDeviceToHost);
  cudaMemcpy(atomForceRec.z, vars->gpu_aForceRecz, sizeof(double) * atomCount, cudaMemcpyDeviceToHost);
  cudaMemcpy(molForceRec.x, vars->gpu_mForceRecx, sizeof(double) * molCount, cudaMemcpyDeviceToHost);
  cudaMemcpy(molForceRec.y, vars->gpu_mForceRecy, sizeof(double) * molCount, cudaMemcpyDeviceToHost);
  cudaMemcpy(molForceRec.z, vars->gpu_mForceRecz, sizeof(double) * molCount, cudaMemcpyDeviceToHost);

  cudaDeviceSynchronize();
  delete [] arr_particleHasNoCharge;
  CUFREE(gpu_particleCharge);
  CUFREE(gpu_particleHasNoCharge);
  CUFREE(gpu_startMol);
  CUFREE(gpu_lengthMol);
  CUFREE(gpu_particleMol);
  CUFREE(gpu_particleKind);
}

__global__ void BoxForceReciprocalGPU(
  double *gpu_aForceRecx,
  double *gpu_aForceRecy,
  double *gpu_aForceRecz,
  double *gpu_mForceRecx,
  double *gpu_mForceRecy,
  double *gpu_mForceRecz,
  double *gpu_particleCharge,
  int *gpu_particleMol,
  int *gpu_particleKind,
  bool *gpu_particleHasNoCharge,
  int *gpu_startMol,
  int *gpu_lengthMol,
  double alpha,
  double alphaSq,
  double qqFact,
  double constValue,
  int imageSize,
  double *gpu_kx,
  double *gpu_ky,
  double *gpu_kz,
  double *gpu_x,
  double *gpu_y,
  double *gpu_z,
  double *gpu_prefact,
  double *gpu_sumRnew,
  double *gpu_sumInew,
  bool *gpu_isFraction,
  int *gpu_molIndex,
  int *gpu_kindIndex,
  double *gpu_lambdaCoulomb,
  double axx,
  double axy,
  double axz,
  int box
)
{
  __shared__ double shared[24];
  int laneID = threadIdx.x % 32;
  int warpID = threadIdx.x / 32;
  int particleID = blockIdx.x;
  double forceX = 0.0, forceY = 0.0, forceZ = 0.0;
  int moleculeID = gpu_particleMol[particleID];
  int kindID = gpu_particleKind[particleID];
  if(!gpu_particleHasNoCharge[particleID]) {
    double lambdaCoef = DeviceGetLambdaCoulomb(moleculeID, kindID, box, gpu_isFraction, gpu_molIndex, gpu_kindIndex, gpu_lambdaCoulomb);

    // loop over other particles within the same molecule
    if(threadIdx.x == 0) {
      double intraForce = 0.0, distSq = 0.0, dist = 0.0;
      double distVectX = 0.0, distVectY = 0.0, distVectZ = 0.0;
      int lastParticleWithinSameMolecule = gpu_startMol[particleID] + gpu_lengthMol[particleID];
      for(int otherParticle = gpu_startMol[particleID];
        otherParticle < lastParticleWithinSameMolecule;
        otherParticle++)
      {
        if(particleID != otherParticle) {
          DeviceInRcut(distSq, distVectX, distVectY, distVectZ, gpu_x, gpu_y, gpu_z, particleID, otherParticle, axx, axy, axz, box);
          dist = sqrt(distSq);

          double expConstValue = exp(-1.0 * alphaSq * distSq);
          double qiqj = gpu_particleCharge[particleID] * gpu_particleCharge[otherParticle] * qqFact;
          intraForce = qiqj * lambdaCoef * lambdaCoef / distSq;
          intraForce *= ((erf(alpha * dist) / dist) - constValue * expConstValue);
          forceX -= intraForce * distVectX;
          forceY -= intraForce * distVectY;
          forceZ -= intraForce * distVectZ;
        }
      }
    }

    // loop over images
    for(int vectorIndex = threadIdx.x; vectorIndex < imageSize; vectorIndex += blockDim.x) {
      double dot = gpu_x[particleID] * gpu_kx[vectorIndex] +
        gpu_y[particleID] * gpu_ky[vectorIndex] + 
        gpu_z[particleID] * gpu_kz[vectorIndex];
        
      double factor = 2.0 * gpu_particleCharge[particleID] * gpu_prefact[vectorIndex] * lambdaCoef *
        (sin(dot) * gpu_sumRnew[vectorIndex] - cos(dot) * gpu_sumInew[vectorIndex]);
        
      forceX += factor * gpu_kx[vectorIndex];
      forceY += factor * gpu_ky[vectorIndex];
      forceZ += factor * gpu_kz[vectorIndex];
    }
  }

  // perform reduction at this point
  int warpSize = 32;
  for (int offset = warpSize/2; offset > 0; offset /= 2) {
    forceX += __shfl_down_sync(FULL_MASK, forceX, offset);
    forceY += __shfl_down_sync(FULL_MASK, forceY, offset);
    forceZ += __shfl_down_sync(FULL_MASK, forceZ, offset);
  }
  if(laneID == 0) {
    shared[warpID*3+0] = forceX;
    shared[warpID*3+1] = forceY;
    shared[warpID*3+2] = forceZ;
  }

  // first thread inside the block will write back to global memory
  __syncthreads();
  if(threadIdx.x == 0) {
    for(int w=1; w<8; w++) {
      forceX += shared[w*3+0];
      forceY += shared[w*3+1];
      forceZ += shared[w*3+2];
    }
    gpu_aForceRecx[particleID] = forceX;
    gpu_aForceRecy[particleID] = forceY;
    gpu_aForceRecz[particleID] = forceZ;
    atomicAdd(&gpu_mForceRecx[moleculeID], forceX);
    atomicAdd(&gpu_mForceRecy[moleculeID], forceY);
    atomicAdd(&gpu_mForceRecz[moleculeID], forceZ);
  }
}

__global__ void SwapReciprocalGPU(double *gpu_x, double *gpu_y, double *gpu_z,
                                  double *gpu_kx, double *gpu_ky, double *gpu_kz,
                                  int atomNumber,
                                  double *gpu_particleCharge,
                                  double *gpu_sumRnew,
                                  double *gpu_sumInew,
                                  double *gpu_sumRref,
                                  double *gpu_sumIref,
                                  double *gpu_prefactRef,
                                  int insert,
                                  double *gpu_energyRecipNew,
                                  int imageSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= imageSize)
    return;

  int p;
  double dotProduct = 0.0, sumReal = 0.0, sumImaginary = 0.0;

  for(p = 0; p < atomNumber; p++) {
    dotProduct = DotProductGPU(gpu_kx[threadID], gpu_ky[threadID],
                               gpu_kz[threadID], gpu_x[p], gpu_y[p], gpu_z[p]);
    sumReal += (gpu_particleCharge[p] * cos(dotProduct));
    sumImaginary += (gpu_particleCharge[p] * sin(dotProduct));
  }

  //If we insert the molecule to the box, we add the sum value.
  //Otherwise, we subtract the sum value
  if(insert) {
    gpu_sumRnew[threadID] = gpu_sumRref[threadID] + sumReal;
    gpu_sumInew[threadID] = gpu_sumIref[threadID] + sumImaginary;
  } else {
    gpu_sumRnew[threadID] = gpu_sumRref[threadID] - sumReal;
    gpu_sumInew[threadID] = gpu_sumIref[threadID] - sumImaginary;
  }

  gpu_energyRecipNew[threadID] = ((gpu_sumRnew[threadID] *
                                   gpu_sumRnew[threadID] +
                                   gpu_sumInew[threadID] *
                                   gpu_sumInew[threadID]) *
                                  gpu_prefactRef[threadID]);
}

__global__ void MolReciprocalGPU(double *gpu_cx, double *gpu_cy, double *gpu_cz,
                                 double *gpu_nx, double *gpu_ny, double *gpu_nz,
                                 double *gpu_kx, double *gpu_ky, double *gpu_kz,
                                 int atomNumber,
                                 double *gpu_particleCharge,
                                 double *gpu_sumRnew,
                                 double *gpu_sumInew,
                                 double *gpu_sumRref,
                                 double *gpu_sumIref,
                                 double *gpu_prefactRef,
                                 double *gpu_energyRecipNew,
                                 int imageSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= imageSize)
    return;
  int p;
  double dotProductOld = 0.0, dotProductNew = 0.0;
  double sumRealNew = 0.0, sumImaginaryNew = 0.0;
  double sumRealOld = 0.0, sumImaginaryOld = 0.0;

  for(p = 0; p < atomNumber; p++) {
    dotProductOld = DotProductGPU(gpu_kx[threadID], gpu_ky[threadID],
                                  gpu_kz[threadID],
                                  gpu_cx[p], gpu_cy[p], gpu_cz[p]);
    dotProductNew = DotProductGPU(gpu_kx[threadID], gpu_ky[threadID],
                                  gpu_kz[threadID],
                                  gpu_nx[p], gpu_ny[p], gpu_nz[p]);
    sumRealNew += (gpu_particleCharge[p] * cos(dotProductNew));
    sumImaginaryNew += (gpu_particleCharge[p] * sin(dotProductNew));
    sumRealOld += (gpu_particleCharge[p] * cos(dotProductOld));
    sumImaginaryOld += (gpu_particleCharge[p] * sin(dotProductOld));
  }

  gpu_sumRnew[threadID] = gpu_sumRref[threadID] - sumRealOld + sumRealNew;
  gpu_sumInew[threadID] = gpu_sumIref[threadID] - sumImaginaryOld +
                          sumImaginaryNew;

  gpu_energyRecipNew[threadID] = ((gpu_sumRnew[threadID] *
                                   gpu_sumRnew[threadID] +
                                   gpu_sumInew[threadID] *
                                   gpu_sumInew[threadID]) *
                                  gpu_prefactRef[threadID]);
}

__global__ void BoxReciprocalSetupGPU(double *gpu_x,
                                      double *gpu_y,
                                      double *gpu_z,
                                      double *gpu_kx,
                                      double *gpu_ky,
                                      double *gpu_kz,
                                      double atomNumber,
                                      double *gpu_particleCharge,
                                      double *gpu_sumRnew,
                                      double *gpu_sumInew,
                                      int imageSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= imageSize)
    return;
  int i;
  double dotP;

  gpu_sumRnew[threadID] = 0.0;
  gpu_sumInew[threadID] = 0.0;
  for(i = 0; i < atomNumber; i++) {
    dotP = DotProductGPU(gpu_kx[threadID], gpu_ky[threadID], gpu_kz[threadID],
                         gpu_x[i], gpu_y[i], gpu_z[i]);
    gpu_sumRnew[threadID] += gpu_particleCharge[i] * cos(dotP);
    gpu_sumInew[threadID] += gpu_particleCharge[i] * sin(dotP);
  }
}

__global__ void BoxReciprocalGPU(double *gpu_prefact,
                                 double *gpu_sumRnew,
                                 double *gpu_sumInew,
                                 double *gpu_energyRecip,
                                 int imageSize)
{
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  if(threadID >= imageSize)
    return;

  gpu_energyRecip[threadID] = ((gpu_sumRnew[threadID] * gpu_sumRnew[threadID] +
                                gpu_sumInew[threadID] * gpu_sumInew[threadID]) *
                               gpu_prefact[threadID]);
}

#endif
