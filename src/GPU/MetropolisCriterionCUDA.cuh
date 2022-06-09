/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at <https://opensource.org/licenses/MIT>.
********************************************************************************/
#pragma once
#ifdef GOMC_CUDA

#include "HelperFunctionsCUDA.cuh"
#include "CalculateMinImageCUDAKernel.cuh"

void CallBMPAccept(VariablesCUDA *vars,
                    int molCount);

__global__ void GetCoeffTranslation(   
                            int numberOfMolecules,
                            double t_max,
                            double * mp_coefficient,
                            double * BETA,
                            double * t_k_x,
                            double * t_k_y,
                            double * t_k_z,
                            double * molForceRefX,
                            double * molForceRefY,
                            double * molForceRefZ,
                            double * molForceNewX,
                            double * molForceNewY,
                            double * molForceNewZ,
                            double * molForceRecRefX,
                            double * molForceRecRefY,
                            double * molForceRecRefZ,
                            double * molForceRecNewX,
                            double * molForceRecNewY,
                            double * molForceRecNewZ);

__device__ double CalculateWRatio(  const double3 &lb_new,
                                    const double3 &lb_old,
                                    const double3 &k,
                                    const double max4);

#endif