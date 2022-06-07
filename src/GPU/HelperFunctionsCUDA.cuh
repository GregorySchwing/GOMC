#pragma once
#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
  
  inline __device__ double LengthSq(double3 & vec){
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
  }
  

#endif