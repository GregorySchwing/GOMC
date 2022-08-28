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

class CUDAMemoryUtils
{
public:
static void ZeroForces(VariablesCUDA *vars,
                     XYZArray const &coords,
                     const MoleculeLookup& molLookup,
                     uint const buffer_index);

};
#endif