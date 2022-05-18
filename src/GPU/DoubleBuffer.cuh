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
//https://stackoverflow.com/questions/2008948/double-buffering-for-game-objects-whats-a-nice-clean-generic-c-way

/**
 * @brief A helper class for accessing parts of a global array to logically
            create multiple buffers.
 * 
 * @tparam T - datatype of array
 */

template< typename T>
class DeviceArray
{
public:
    DeviceArray(){};

    T *  get() const { return _n; }
    void set(T * n) { _n = n; }

private:
    T * _n;
};

/**
 * @brief   Class for managing a general number of states of a single array.
            Allocates one array and logically partitions it into multiple arrays.
            Also, manages the buffer states.
 * 
 * @tparam T - helper class to hold partial array
 * @tparam DAT - datatype
 * @tparam n - number of states
 */
template< class T, typename DAT, std::size_t n >
struct MultiBuffer
{
    MultiBuffer() : _active_offset(0) {}
    MultiBuffer(std::size_t m) : _active_offset(0) {
        CUMALLOC((void**) &_globalmemory, m * sizeof(DAT));
        for (int index = 0; index < n; ++index)
            _objects[index].set(&_globalmemory[index*m]);
    }

    void ChangeBuffers() { ++_active_offset; }
    T* GetInstance(std::size_t k) { return &_objects[ (_active_offset + k) % n ]; }

private:
    T _objects[n];
    std::size_t _active_offset;
    DAT * _globalmemory;
};

template< class T, typename DAT, std::size_t n >
class BufferAccess
{
public:
    BufferAccess( MultiBuffer< T, DAT, n >& buf, std::size_t offset )
        : _buffer(buf), _offset(offset)
    {
    }

    T* operator->() const
    {
        return _buffer.GetInstance(_offset);
    }

private:
    MultiBuffer< T, DAT, n >& _buffer;
    const std::size_t _offset;
};
#endif
#endif