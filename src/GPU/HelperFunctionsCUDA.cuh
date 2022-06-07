#pragma once
#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
  
    inline __device__ double LengthSq(double3 & vec){
        return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
    }

    inline __device__ double3 Subtract(double3 & vec1, double3 & vec2){
        return make_double3(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z);
    }
  
    inline __device__ double3 Add(double3 & vec1, double3 & vec2){
        return make_double3(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z);
    }

#endif