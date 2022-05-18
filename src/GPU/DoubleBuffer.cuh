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
/*
template <typename Buf_Type>
class DoubleBuffer {
        DoubleBuffer();
        Buf_Type ** selector;
        int size;
        Buf_Type * array1;
        Buf_Type * array2;
    public:
        DoubleBuffer<Buf_Type>(int);
        ~DoubleBuffer();
};
*/
template< class T>
class TestClass
{
public:
    TestClass() : _n(0) {}

    int get() const { return _n; }
    void set(int n) { _n = n; }

private:
    T * _n;
};
template< class T, std::size_t n >
struct MultiBuffer
{
    MultiBuffer() : _active_offset(0) {}
    MultiBuffer(std::size_t m) : _active_offset(0) {}

    void ChangeBuffers() { ++_active_offset; }
    T* GetInstance(std::size_t k) { return &_objects[ (_active_offset + k) % n ]; }

private:
    T _objects[n];
    std::size_t _active_offset;
};

template< class T, std::size_t n >
class BufferAccess
{
public:
    BufferAccess( MultiBuffer< T, n >& buf, std::size_t offset )
        : _buffer(buf), _offset(offset)
    {
    }

    T* operator->() const
    {
        return _buffer.GetInstance(_offset);
    }

private:
    MultiBuffer< T, n >& _buffer;
    const std::size_t _offset;
};
#include "DoubleBuffer.cu"
#endif
#endif