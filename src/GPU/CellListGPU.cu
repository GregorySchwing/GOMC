#ifdef GOMC_CUDA
#include "CellListGPU.cuh"
#include "cub/cub.cuh"

CellListGPU::CellListGPU(VariablesCUDA * cv, int _atomNumber, MoleculeLookup & molLookup) : 
atomNumber(_atomNumber), molLookRef(molLookup)
{
    CUMALLOC((void**) &cv->gpu_Ones, atomNumber * sizeof(int));
    CUMALLOC((void**) &cv->gpu_particleIndices, atomNumber * sizeof(int));
    pI.resize(atomNumber);
	thrust::sequence(pI.begin(),pI.end());
    ones.resize(atomNumber);
    thrust::fill(ones.begin(),ones.end(), 1);
	thrust::sequence(pI.begin(),pI.end());
    // Fill with 1s
    cudaMemcpy(cv->gpu_Ones, thrust::raw_pointer_cast(&ones[0]), atomNumber * sizeof(int), cudaMemcpyDeviceToDevice);
    cudaMemcpy(cv->gpu_particleIndices, thrust::raw_pointer_cast(&pI[0]), atomNumber * sizeof(int), cudaMemcpyDeviceToDevice);
}
/*
void CellListGPU::GridBox(VariablesCUDA * cv,
                        XYZArray const &coords,
                        XYZArray const &axes,
                        const int buffer_index,
                        const uint b){
    GOMC_EVENT_START(1, GomcProfileEvent::GRID_ALL_GPU);

    int atomCount = molLookRef.NumAtomsInBox(b);
    int boxAtomOffset = b ? atomCount : 0;

    // Need to reinitialize the sequence 0..N-1 every GridAll
    cudaMemcpy(cv->gpu_particleIndices, thrust::raw_pointer_cast(&pI[0]), atomNumber * sizeof(int), cudaMemcpyDeviceToDevice);
    // Clear Cell Degrees
    //cuMemsetD32(reinterpret_cast<CUdeviceptr>(cv->gpu_cellDegrees),  0, size_t(numberOfCells));
    

    // Make sure cpu_numberOfCells == gpu_numberOfCells
    cudaMemset(&cv->gpu_cellDegrees[cv->cpu_startOfBoxCellList[b]], 0, cv->cpu_numberOfCells[b]*sizeof(int));

    BufferAccess<DeviceArray<int>, int, buffers> mapParticleToCell_view(*(cv->gpu_mapParticleToCell), buffer_index);
    BufferAccess<DeviceArray<int>, int, buffers> cellVector_view(*(cv->gpu_cellVector), buffer_index);
    BufferAccess<DeviceArray<int>, int, buffers> cellStartIndex_view(*(cv->gpu_cellStartIndex), buffer_index);

    BufferAccess<DeviceArray<double>, double, buffers> coords_x_view(*(cv->gpu_coords_x), buffer_index);
    BufferAccess<DeviceArray<double>, double, buffers> coords_y_view(*(cv->gpu_coords_y), buffer_index);
    BufferAccess<DeviceArray<double>, double, buffers> coords_z_view(*(cv->gpu_coords_z), buffer_index);

    MapParticlesToCell(cv,
                    coords_x_view->get(),
                    coords_y_view->get(),
                    coords_z_view->get(),  
                    mapParticleToCell_view->get(),  
                    atomCount,
                    axes,
                    b);
    SortMappedParticles(cv,
                        mapParticleToCell_view->get(),  
                        cellVector_view->get(),  
                        coords);

    CalculateCellDegrees(cv,coords);
    PrefixScanCellDegrees(cv, cellStartIndex_view->get(), cv->cpu_numberOfCells[b]);
    GOMC_EVENT_STOP(1, GomcProfileEvent::GRID_ALL_GPU);
}
*/
void CellListGPU::GridAll(VariablesCUDA * cv,
                        XYZArray const &coords,
                        XYZArray const &axes,
                        int numberOfCells,
                        const int buffer_index){
    GOMC_EVENT_START(1, GomcProfileEvent::GRID_ALL_GPU);
    printf("GridAll CLGPU\n");

    int atomCount = coords.Count();

    // Need to reinitialize the sequence 0..N-1 every GridAll
    cudaMemcpy(cv->gpu_particleIndices, thrust::raw_pointer_cast(&pI[0]), atomNumber * sizeof(int), cudaMemcpyDeviceToDevice);
    // Clear Cell Degrees
    //cuMemsetD32(reinterpret_cast<CUdeviceptr>(cv->gpu_cellDegrees),  0, size_t(numberOfCells));
    cudaMemset(cv->gpu_cellDegrees, 0, numberOfCells*sizeof(int));


    BufferAccess<DeviceArray<int>, int, buffers> mapParticleToCell_view(*(cv->gpu_mapParticleToCell), buffer_index);
    BufferAccess<DeviceArray<int>, int, buffers> cellVector_view(*(cv->gpu_cellVector), buffer_index);
    BufferAccess<DeviceArray<int>, int, buffers> cellStartIndex_view(*(cv->gpu_cellStartIndex), buffer_index);

    BufferAccess<DeviceArray<double>, double, buffers> coords_x_view(*(cv->gpu_coords_x), buffer_index);
    BufferAccess<DeviceArray<double>, double, buffers> coords_y_view(*(cv->gpu_coords_y), buffer_index);
    BufferAccess<DeviceArray<double>, double, buffers> coords_z_view(*(cv->gpu_coords_z), buffer_index);
    MapParticlesToCell(cv,
                    coords_x_view->get(),
                    coords_y_view->get(),
                    coords_z_view->get(),  
                    mapParticleToCell_view->get(),  
                    axes,
                    molLookRef);
    SortMappedParticles(cv,
                        mapParticleToCell_view->get(),  
                        cellVector_view->get(),  
                        coords);
    CalculateCellDegrees(cv,coords);
    PrefixScanCellDegrees(cv, cellStartIndex_view->get(), numberOfCells);
    GOMC_EVENT_STOP(1, GomcProfileEvent::GRID_ALL_GPU);
}


void CellListGPU::MapParticlesToCell(VariablesCUDA * cv,
                                    double * x,
                                    double * y,
                                    double * z,
                                    int * mp2c,
                                    XYZArray const &axes,
                                    MoleculeLookup & molLookupRef){
    // Run the kernel
    int threadsPerBlock = 256;
    int molCount = molLookupRef.molLookupCount;
    int blocksPerGrid = (int)((molCount * warp_size) / threadsPerBlock) + 1;

    MapParticlesToCellKernel<<< blocksPerGrid, threadsPerBlock>>>(
                            molCount,
                            x,
                            y,
                            z,                               
                            mp2c,
                            molLookupRef.molLookupGPU->gpu_molLookup,
                            molLookupRef.molLookupGPU->gpu_numMolsInBox,
                            molLookupRef.molLookupGPU->gpu_startAtomIdx,
                            cv->gpu_cellSize,
                            cv->gpu_edgeCells,
                            cv->gpu_nonOrth,
                            cv->gpu_Invcell_x,
                            cv->gpu_Invcell_y,
                            cv->gpu_Invcell_z);
    cudaDeviceSynchronize();
    checkLastErrorCUDA(__FILE__, __LINE__);

}


void CellListGPU::SortMappedParticles(VariablesCUDA * cv,
                                    int * mp2c,
                                    int * cellVec,
                                    XYZArray const &coords){
    int atomNumber = coords.Count();
    
    // Run the kernel
    CreateStartVector(atomNumber,
                    mp2c,
                    cv->gpu_mapParticleToCellSorted,
                    cv->gpu_particleIndices,
                    cellVec,
                    cv->d_temp_storage_sort,
                    cv->temp_storage_bytes_sort);
    cudaDeviceSynchronize();
    checkLastErrorCUDA(__FILE__, __LINE__);

}

void CellListGPU::CalculateCellDegrees(VariablesCUDA * cv,
                                    XYZArray const &coords){
    int atomNumber = coords.Count();                                            // Run the kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (int)(atomNumber / threadsPerBlock) + 1;
    // Run the kernel
    /*
    CalculateCellDegreesKernel<<< blocksPerGrid, 
                                threadsPerBlock, 
                                (3*threadsPerBlock+2)*sizeof(int)>>>(atomNumber,
                                                                cv->gpu_mapParticleToCellSorted,
                                                                cv->gpu_cellDegrees);
    */
 CalculateCellDegreesKernel<<< blocksPerGrid, 
                                threadsPerBlock>>>(atomNumber,
                                                    cv->gpu_mapParticleToCellSorted,
                                                    cv->gpu_cellDegrees);
    cudaDeviceSynchronize();
    checkLastErrorCUDA(__FILE__, __LINE__);

}

void CellListGPU::CalculateCellDegreesCUB(VariablesCUDA * cv,
                                    XYZArray const &coords){
    int atomNumber = coords.Count();                                            // Run the kernel
    CreateCellDegrees(atomNumber,
                    cv->gpu_mapParticleToCellSorted,
                    cv->gpu_Ones,
                    cv->gpu_CellDegreeSanityCheck,
                    cv->gpu_cellDegrees,
                    cv->gpu_IterationsReq);
    cudaDeviceSynchronize();
    checkLastErrorCUDA(__FILE__, __LINE__);

}
// CustomMin functor
struct CustomMin
{
    template <typename T>
    CUB_RUNTIME_FUNCTION __forceinline__
    T operator()(const T &a, const T &b) const {
        return (b < a) ? b : a;
    }
};
// Need to sort first, since ReduceByKey needs keys to be in stretches of the same key
void CellListGPU::CreateCellDegrees(int numberOfAtoms,
                                int * mapParticleToCellSorted,
                                int * Ones,
                                int * gpu_CellDegreeSanityCheck,
                                int * cellDegrees,
                                int * iterationsReq){
    // Declare, allocate, and initialize device-accessible pointers for input and output
    int          num_items = numberOfAtoms;          // e.g., 8
    int          *d_keys_in = mapParticleToCellSorted;           // e.g., [0, 2, 2, 9, 5, 5, 5, 8]
    int          *d_values_in = Ones;                            // e.g., [0, 7, 1, 6, 2, 5, 3, 4]
    int          *d_unique_out = gpu_CellDegreeSanityCheck;      // e.g., [-, -, -, -, -, -, -, -]
    int          *d_aggregates_out = cellDegrees;  // e.g., [-, -, -, -, -, -, -, -]
    int          *d_num_runs_out = iterationsReq;    // e.g., [-]
    CustomMin    reduction_op;
    // Determine temporary device storage requirements
    void     *d_temp_storage = NULL;
    size_t   temp_storage_bytes = 0;
    cub::DeviceReduce::ReduceByKey(d_temp_storage, temp_storage_bytes, d_keys_in, d_unique_out, d_values_in, d_aggregates_out, d_num_runs_out, reduction_op, num_items);
    // Allocate temporary storage
    cudaMalloc(&d_temp_storage, temp_storage_bytes);
    // Run reduce-by-key
    cub::DeviceReduce::ReduceByKey(d_temp_storage, temp_storage_bytes, d_keys_in, d_unique_out, d_values_in, d_aggregates_out, d_num_runs_out, reduction_op, num_items);
    // d_unique_out      <-- [0, 2, 9, 5, 8]
    // d_aggregates_out  <-- [0, 1, 6, 2, 4]
    // d_num_runs_out    <-- [5]
}


void CellListGPU::PrefixScanCellDegrees(VariablesCUDA * cv,
                                    int *csi,
                                    int numberOfCells){
    CalculateNewRowOffsets(numberOfCells,
                           csi,
                           cv->gpu_cellDegrees,
                           cv->d_temp_storage_prefix_sum,
                           cv->temp_storage_bytes_prefix_sum);
    cudaDeviceSynchronize();
    checkLastErrorCUDA(__FILE__, __LINE__);
}

void CellListGPU::CopyGPUMemoryToToHost(int * deviceMemory,
                                    int size,
                                    std::vector<int> & hostMemory){
    hostMemory.clear();
    hostMemory.resize(size);
    cudaMemcpy(&hostMemory[0], deviceMemory, size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    checkLastErrorCUDA(__FILE__, __LINE__);
}

__device__ int PositionToCell(int atomIndex,
                            double* gpu_x,
                            double* gpu_y,
                            double* gpu_z,                                
                            double* gpu_cellSize,
                            int* gpu_edgeCells,
                            int* gpu_nonOrth,
                            double *gpu_Invcell_x,
                            double *gpu_Invcell_y,
                            double *gpu_Invcell_z,
                            const int b = 0){
    double3 pos = make_double3(gpu_x[atomIndex],
                               gpu_y[atomIndex],
                               gpu_z[atomIndex]);
    if(gpu_nonOrth[0]){
        pos = TransformUnSlantGPU(pos, 
                                gpu_Invcell_x,
                                gpu_Invcell_y,
                                gpu_Invcell_z);
    }
    int x = (int)(pos.x / gpu_cellSize[3*b + 0]);
    int y = (int)(pos.y / gpu_cellSize[3*b + 1]);
    int z = (int)(pos.z / gpu_cellSize[3*b + 2]);
    //Check the cell number to avoid segfaults for coordinates close to axis
    //x, y, and z should never be equal or greater than number of cells in x, y,
    // and z axis, respectively.
    x -= (x == gpu_edgeCells[3*b + 0] ?  1 : 0);
    y -= (y == gpu_edgeCells[3*b + 1] ?  1 : 0);
    z -= (z == gpu_edgeCells[3*b + 2] ?  1 : 0);
    return x * gpu_edgeCells[1] * gpu_edgeCells[2] + y * gpu_edgeCells[2] + z;
}

__global__ void MapParticlesToCellKernel(
                            int molCount,
                            double* gpu_x,
                            double* gpu_y,
                            double* gpu_z,                                
                            int* gpu_mapParticleToCell,
                            uint* gpu_molLookup,
                            uint* gpu_molBoxCount,
                            int* gpu_startAtomIdx,
                            double *gpu_cellSize,
                            int *gpu_edgeCells,
                            int* gpu_nonOrth,
                            double **gpu_Invcell_x,
                            double **gpu_Invcell_y,
                            double **gpu_Invcell_z){
    int threadID = blockIdx.x * blockDim.x + threadIdx.x;
    // Optimal alu usage/memory latency will probably be 1 warp/molecule
    uint warpIdx = threadID / warp_size;
    uint laneIdx = threadID % warp_size;
    if (warpIdx >= molCount)
        return;

    int molIndex = gpu_molLookup[warpIdx];
    uint b = molIndex < gpu_molBoxCount[0];
    //printf("b %d", b);
    for (int particleIndex = gpu_startAtomIdx[molIndex] + laneIdx; particleIndex < gpu_startAtomIdx[molIndex + 1]; particleIndex += warp_size ){
        int cell = PositionToCell(particleIndex,
                                gpu_x,
                                gpu_y,
                                gpu_z,
                                gpu_cellSize,
                                gpu_edgeCells,
                                gpu_nonOrth,
                                gpu_Invcell_x[b],
                                gpu_Invcell_y[b],
                                gpu_Invcell_z[b]);
        gpu_mapParticleToCell[particleIndex] = cell;
    }
}



// Make sure a zero element is padded onto the end.
// https://github.com/NVIDIA/cub/issues/367
void CellListGPU::CalculateNewRowOffsets( int numberOfRows,
                                        int * global_row_offsets_dev_ptr,
                                        int * global_degrees_dev_ptr,
                                        void     **d_temp_storage_prefix_sum,
                                        size_t   *temp_storage_bytes_prefix_sum){
    // Declare, allocate, and initialize device-accessible pointers for input and output
    int  num_items = numberOfRows+1;      // e.g., 7
    int  *d_in = global_degrees_dev_ptr;        // e.g., [8, 6, 7, 5, 3, 0, 9]
    int  *d_out = global_row_offsets_dev_ptr;         // e.g., [ ,  ,  ,  ,  ,  ,  ]
    // Determine temporary device storage requirements
    if(*temp_storage_bytes_prefix_sum == 0){

        cub::DeviceScan::ExclusiveSum(*d_temp_storage_prefix_sum, *temp_storage_bytes_prefix_sum, d_in, d_out, num_items);
        // Allocate temporary storage
        CUMALLOC(d_temp_storage_prefix_sum, *temp_storage_bytes_prefix_sum);
    } else {
        cudaMemset(*d_temp_storage_prefix_sum, 0, *temp_storage_bytes_prefix_sum);
    }
    // Run exclusive prefix sum
    cub::DeviceScan::ExclusiveSum(*d_temp_storage_prefix_sum, *temp_storage_bytes_prefix_sum, d_in, d_out, num_items);
    // d_out s<-- [0, 8, 14, 21, 26, 29, 29]
    //cudaFree(temp_storage_bytes_prefix_sum);
}



void CellListGPU::CreateStartVector(int numberOfAtoms,
                                int * mapParticleToCell,
                                int * mapParticleToCellSorted,
                                int * particleIndices,
                                int * particleIndicesSorted,
                                void     **d_temp_storage_sort,
                                size_t   *temp_storage_bytes_sort){
    // Declare, allocate, and initialize device-accessible pointers for sorting data
    // numberOfAtoms            e.g., 7
    // mapParticleToCell        e.g., [8, 6, 7, 5, 3, 0, 9]
    // mapParticleToCellSorted  e.g., [        ...        ]
    // particleIndices          e.g., [0, 1, 2, 3, 4, 5, 6]
    // particleIndicesSorted    e.g., [        ...        ]
    // Determine temporary device storage requirements
    int num_items = numberOfAtoms;
    int  *d_keys_in = mapParticleToCell;
    int  *d_keys_out = mapParticleToCellSorted;       
    int  *d_values_in = particleIndices;    
    int  *d_values_out = particleIndicesSorted;   
    if(*temp_storage_bytes_sort == 0){
        cub::DeviceRadixSort::SortPairs(*d_temp_storage_sort, *temp_storage_bytes_sort,
            d_keys_in, d_keys_out, d_values_in, d_values_out, num_items);
        // Allocate temporary storage
        CUMALLOC(d_temp_storage_sort, *temp_storage_bytes_sort);
    } else {
        cudaMemset(*d_temp_storage_sort, 0, *temp_storage_bytes_sort);
    }
    // Run sorting operation
    cub::DeviceRadixSort::SortPairs(*d_temp_storage_sort, *temp_storage_bytes_sort,
        d_keys_in, d_keys_out, d_values_in, d_values_out, num_items);
    // mapParticleToCellSorted        <-- [0, 3, 5, 6, 7, 8, 9]
    // particleIndicesSorted          <-- [5, 4, 3, 1, 2, 0, 6]
}

__global__ void CalculateCellDegreesKernel(int atomNumber,
                                            int* gpu_mapParticleToCellSorted,
                                            int* gpu_cellDegrees){
    int threadID = blockIdx.x * blockDim.x + threadIdx.x;
    if (threadID >= atomNumber)
        return;
    atomicAdd(&gpu_cellDegrees[gpu_mapParticleToCellSorted[threadID]], 1);
}

/*
__global__ void CalculateCellDegreesKernel(int atomNumber,
                                            int* gpu_mapParticleToCellSorted,
                                            int* gpu_cellDegrees){
    int threadID = blockIdx.x * blockDim.x + threadIdx.x;
    if (threadID >= atomNumber)
        return;
    extern __shared__ int part_ary[];
    // minimum cell to maximum cell
    // 0..BlockSize-1 - count of each cell key
    // BlockSize..2*BlockSize-1 - indicator variable (0,1)
    // 2*BlockSize..3*BlockSize-1 - load cell key from global memory
    // 3*BlockSize - Min cell key
    // 3*BlockSize+1 - Max cell key
    part_ary[2*blockDim.x+threadIdx.x] = gpu_mapParticleToCellSorted[threadID];
    part_ary[threadIdx.x] = part_ary[2*blockDim.x+threadIdx.x];
    int i = blockDim.x/2;
    while (0 < i){
        if(threadIdx.x < i){
            part_ary[threadIdx.x] = min(part_ary[threadIdx.x],part_ary[threadIdx.x+i]); 
        }
          __syncthreads();
        i /= 2;
    }
    if (threadIdx.x==0)
        part_ary[3*blockDim.x] = part_ary[threadIdx.x];

      __syncthreads();

    part_ary[threadIdx.x] = part_ary[2*blockDim.x+threadIdx.x];
    i = blockDim.x/2;
    while (0 < i){
        if(threadIdx.x < i){
            part_ary[threadIdx.x] = max(part_ary[threadIdx.x],part_ary[threadIdx.x+i]); 
        }
          __syncthreads();
        i /= 2;
    }
    if (threadIdx.x==0)
        part_ary[3*blockDim.x+1] = part_ary[threadIdx.x];

      __syncthreads();

    for (int reductionIteration = 0; reductionIteration < part_ary[3*blockDim.x+1]-part_ary[3*blockDim.x]; ++reductionIteration){
        part_ary[blockDim.x+threadIdx.x] =  part_ary[2*blockDim.x+threadIdx.x] == part_ary[3*blockDim.x]+reductionIteration;
        i = blockDim.x/2;
        while (0 < i){
            if(threadIdx.x < i){
                part_ary[blockDim.x+threadIdx.x] = part_ary[blockDim.x+threadIdx.x] + part_ary[blockDim.x+threadIdx.x+i]; 
            }
              __syncthreads();
            i /= 2;
        }
        if (threadIdx.x==0)
            part_ary[reductionIteration] = part_ary[blockDim.x+threadIdx.x];

          __syncthreads();
    }

    for (int reductionIteration = threadIdx.x; reductionIteration < part_ary[3*blockDim.x+1]-part_ary[3*blockDim.x]; ++reductionIteration){
        atomicAdd(&gpu_cellDegrees[part_ary[3*blockDim.x]+reductionIteration], part_ary[reductionIteration]);
    }
}
*/

#endif
