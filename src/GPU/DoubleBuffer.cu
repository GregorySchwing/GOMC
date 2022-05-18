/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at <https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef DOUBLE_BUFFER_CUH
#define DOUBLE_BUFFER_CUH
#include "DoubleBuffer.cuh"
#ifdef GOMC_CUDA

template<typename Buf_Type> 
DoubleBuffer<Buf_Type>::DoubleBuffer<Buf_Type>(int arraySize) {
    CUMALLOC((void**) selector, 2 * sizeof(void*));
    CUMALLOC((void**) array1, arraySize * sizeof(Buf_Type));
    CUMALLOC((void**) array2, arraySize * sizeof(Buf_Type));
    
    cudaMemcpy(&selector[0], array1, 1 * sizeof(void*),
    cudaMemcpyDeviceToDevice);
    cudaMemcpy(&selector[1], array2, 1 * sizeof(void*),
    cudaMemcpyDeviceToDevice);
}
template<typename Buf_Type> DoubleBuffer<Buf_Type>::~DoubleBuffer<Buf_Type>() {
    CUFREE(array1);
    CUFREE(array2);
}
  

#endif
#endif