#ifndef CELLLIST_GPU_H
#define CELLLIST_GPU_H
#ifdef GOMC_CUDA

#include "VariablesCUDA.cuh"
#include "MoleculeLookup.h"
#include "XYZArray.h"
#include "CalculateMinImageCUDAKernel.cuh"
#include<thrust/device_vector.h>
#include<thrust/sequence.h>
#include<thrust/fill.h>
#include "CUDAMemoryManager.cuh"
#include "GOMCEventsProfile.h" // for NVTX profiling


static const int warp_size = 32; 
// For cuMemsetD32
//#include <cuda_runtime.h>

__global__ void MapParticlesToCellKernel(int molCount,
                    double* gpu_x,
                    double* gpu_y,
                    double* gpu_z,                                
                    int* gpu_mapParticleToCell,
                    uint* gpu_molLookup,
                    uint * cv_gpu_molLookup,
                    uint* gpu_molBoxCount,
                    int* gpu_startAtomIdx,                    
                    double *gpu_cellSize,
                    int *gpu_edgeCells,
                    int* gpu_nonOrth,
                    double **gpu_Invcell_x,
                    double **gpu_Invcell_y,
                    double **gpu_Invcell_z);

__global__ void CalculateCellDegreesKernel(int atomNumber,
                                            int* gpu_mapParticleToCellSorted,
                                            int* gpu_cellDegrees);

class CellListGPU {
  public:
    CellListGPU(VariablesCUDA * cv, int atomCount, MoleculeLookup & molLookup);
/*
    void GridBox(VariablesCUDA * cv,
                        XYZArray const &coords,
                        XYZArray const &axes,
                        const int buffer_index,
                        const uint b);
*/
    void GridAll(VariablesCUDA * cv,
                  XYZArray const &coords,
                  XYZArray const &axes,
                  int numberOfCells,
                  const int buffer_index = 0);
    void CopyGPUMemoryToToHost(int * deviceMemory,
                                    int size,
                                    std::vector<int> & hostMemory);

    void MapParticlesToCell(VariablesCUDA * cv,
                                    double * x,
                                    double * y,
                                    double * z,
                                    int * mp2c,
                                    XYZArray const &axes,
                                    MoleculeLookup & molLookupRef);

    void SortMappedParticles(VariablesCUDA * cv,
                                    int * mp2c,
                                    int * cellVec,
                                    XYZArray const &coords);
    void CalculateCellDegrees(VariablesCUDA * cv,
                              XYZArray const &coords);
    void CalculateCellDegreesCUB(VariablesCUDA * cv,
                              XYZArray const &coords);
    void PrefixScanCellDegrees(VariablesCUDA * cv,
                                        int *csi,
                                        int numberOfCells);
  private:
    int atomNumber;
    MoleculeLookup & molLookRef;
    thrust::device_vector<int> pI;
    thrust::device_vector<int> ones;
    void CreateStartVector(int numberOfAtoms,
                          int * mapParticleToCell,
                          int * mapParticleToCellSorted,
                          int * particleIndices,
                          int * particleIndicesSorted,
                          void     **d_temp_storage_sort,
                          size_t   *temp_storage_bytes_sort);

    void CreateCellDegrees(int numberOfAtoms,
                          int * mapParticleToCellSortedGPURes,
                          int * OnesGPURes,
                          int * gpu_CellDegreeSanityCheck,
                          int * cellDegreesGPURes,
                          int * iterationsReq);

    void CalculateNewRowOffsets( int numberOfRows,
                                int * global_row_offsets_dev_ptr,
                                int * global_degrees_dev_ptr,
                                void     **d_temp_storage_prefix_sum,
                                size_t   *temp_storage_bytes_prefix_sum);
};

#endif
#endif