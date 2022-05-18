/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at <https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef DOUBLE_BUFFER_CUH
#define DOUBLE_BUFFER_CUH

#ifdef GOMC_CUDA
#include "CUDAMemoryManager.cuh"

template <typename Buf_Type>
class DoubleBuffer {
        DoubleBuffer();
        Buf_Type ** selector;
        int size;
        Buf_Type * array1;
        Buf_Type * array2;
    public:
        DoubleBuffer(int);
        ~DoubleBuffer();
};

#include "DoubleBuffer.cu"
#endif
#endif