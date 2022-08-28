/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at <https://opensource.org/licenses/MIT>.
********************************************************************************/
#pragma once
#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <unordered_map>
#include <iostream>
#include "VariablesCUDA.cuh"
#include "MoleculeLookup.h"
#include "XYZArray.h"

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
                            double *gpu_mTorquez);

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
                            double *gpu_mTorquez_new);


class CUDAMemoryUtils
{
public:

static void ZeroForces(VariablesCUDA *vars,
                     XYZArray const &coords,
                     const MoleculeLookup& molLookup,
                     uint const buffer_index);
static void CallZeroBoxForces(VariablesCUDA *vars,
                     XYZArray const &coords,
                     const MoleculeLookup& molLookup,
                     uint const buffer_index,
                     uint const box);
static void CallCopyBoxForces(VariablesCUDA *vars,
                     XYZArray const &coords,
                     const MoleculeLookup& molLookup,
                     uint const buffer_index,
                     uint const box);
    private:
    static const int currentStateBufferIndex = 0;
    static const int nextStateBufferIndex = 1;
};
#endif