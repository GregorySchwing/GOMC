#pragma once
#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
  
  inline double LengthSq(double3 & vec) const
  {
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
  }

#endif