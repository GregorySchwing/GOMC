/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at <https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "DoubleBuffer.cuh"
#ifdef GOMC_CUDA

template<typename T> DoubleBuffer<T>::DoubleBuffer(int arraySize) {
    CUMALLOC((void**) array1, arraySize * sizeof(T));
    CUMALLOC((void**) array2, arraySize * sizeof(T));
}
template<typename T> DoubleBuffer<T>::~DoubleBuffer() {
    CUFREE(array1);
    CUFREE(array2);
}
  

#endif