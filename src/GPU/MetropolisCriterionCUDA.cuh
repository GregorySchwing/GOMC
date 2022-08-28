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
#include "MoleculeLookup.h"

const int MPDISPLACE = 0;
const int MPROTATE = 1;

void CallGetCoeff(VariablesCUDA *vars,
                    int moveType,
                    int molCount,
                    double * MPCoeff);

void CallGetCoeff(VariablesCUDA *vars,
                    int moveType,
                    int box,
                    const MoleculeLookup& molLookup,
                    double * MPCoeff);

__global__ void GetCoeffTranslation(   
                            int numberOfMolecules,
                            double * t_max,
                            double * BETA,
                            double * mp_coefficient,
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


__global__ void GetCoeffTranslationBox(   
                            int numberOfMolecules,
                            double * t_max,
                            double * BETA,
                            double * mp_coefficient,
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

__global__ void GetCoeffRotationBox(   
                            int numberOfMolecules,
                            double * r_max,
                            double * BETA,
                            double * mp_coefficient,
                            double * r_k_x,
                            double * r_k_y,
                            double * r_k_z,
                            double * molTorqueRefX,
                            double * molTorqueRefY,
                            double * molTorqueRefZ,
                            double * molTorqueNewX,
                            double * molTorqueNewY,
                            double * molTorqueNewZ);

__global__ void GetCoeffRotation(   
                            int numberOfMolecules,
                            double * r_max,
                            double * BETA,
                            double * mp_coefficient,
                            double * r_k_x,
                            double * r_k_y,
                            double * r_k_z,
                            double * molTorqueRefX,
                            double * molTorqueRefY,
                            double * molTorqueRefZ,
                            double * molTorqueNewX,
                            double * molTorqueNewY,
                            double * molTorqueNewZ);
/*
__global__ void Accept(   
                            double mp_coefficient,
                            double BETA,
                            double sysPotNew,
                            double sysPotRef);
*/
__device__ double CalculateWRatio(  const double3 &lb_new,
                                    const double3 &lb_old,
                                    const double3 &k,
                                    const double max4);

#endif