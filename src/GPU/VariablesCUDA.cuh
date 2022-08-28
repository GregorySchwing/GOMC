/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at <https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifdef GOMC_CUDA
#ifndef VARIABLES_CUDA_CUH
#define VARIABLES_CUDA_CUH
#include <cuda.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include "EnsemblePreprocessor.h"
#include "NumLib.h"
#include "DoubleBuffer.cuh"

static const int warp_size = 32; 


//Need a separate float constant for device code with the MSVC compiler
//See CUDA Programming Guide section I.4.13 for details 
static const __device__ double qqFactGPU = num::qqFact;
// Number of buffers in a multi-buffer.  For now only 2.
// Could eventually use one force to try > 1 sequential coord states
const std::size_t buffers = 2;
const int currentState = 0;
const int nextState = 1;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

inline void checkLastErrorCUDA(const char *file, int line)
{
  cudaError_t code = cudaGetLastError();
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    exit(code);
  }
}

inline void printFreeMemory()
{
  size_t free_byte ;
  size_t total_byte ;
  cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;

  if ( cudaSuccess != cuda_status ) {
    printf("Error: cudaMemGetInfo fails, %s \n",
           cudaGetErrorString(cuda_status) );
    exit(1);
  }
  double free_db = (double)free_byte ;
  double total_db = (double)total_byte ;
  double used_db = total_db - free_db ;
  printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
         used_db / 1024.0 / 1024.0, free_db / 1024.0 / 1024.0, total_db / 1024.0 / 1024.0);
}

class VariablesCUDA
{
public:
  VariablesCUDA()
  {
    gpu_sigmaSq = NULL;
    gpu_epsilon_Cn = NULL;
    gpu_n = NULL;
    gpu_VDW_Kind = NULL;
    gpu_isMartini = NULL;
    gpu_count = NULL;
    gpu_rCut = NULL;
    gpu_rCutLow = NULL;
    gpu_rOn = NULL;
    gpu_alpha = NULL;
    gpu_rCutCoulomb = NULL;
    gpu_ewald = NULL;
    gpu_diElectric_1 = NULL;
    
    gpu_aFx = NULL;
    gpu_aFy = NULL;
    gpu_aFz = NULL;
    gpu_mFx = NULL;
    gpu_mFy = NULL;
    gpu_mFz = NULL;

    gpu_aForcex = NULL;
    gpu_aForcey = NULL;
    gpu_aForcez = NULL;
    gpu_mForcex = NULL;
    gpu_mForcey = NULL;
    gpu_mForcez = NULL;

    gpu_aForcex_buffer = NULL;
    gpu_aForcey_buffer = NULL;
    gpu_aForcez_buffer = NULL;
    gpu_mForcex_buffer = NULL;
    gpu_mForcey_buffer = NULL;
    gpu_mForcez_buffer = NULL;
    gpu_startAtomIdx = NULL;

    // setting lambda values to null
    gpu_molIndex = NULL;
    gpu_lambdaVDW = NULL;
    gpu_lambdaCoulomb = NULL;
    gpu_isFraction = NULL;
  }
  double *gpu_sigmaSq;
  double *gpu_epsilon_Cn;
  double *gpu_n;
  int *gpu_VDW_Kind;
  int *gpu_isMartini;
  int *gpu_count;
  int *gpu_startAtomIdx; //start atom index of the molecule
  double *gpu_rCut;
  double *gpu_rCutCoulomb;
  double *gpu_rCutLow;
  double *gpu_rOn;
  double *gpu_alpha;
  int *gpu_ewald;
  double *gpu_diElectric_1;
  double *gpu_x, *gpu_y, *gpu_z;
  // Single molecule arrays
  double *gpu_cx, *gpu_cy, *gpu_cz;
  double *gpu_nx, *gpu_ny, *gpu_nz;
  double *gpu_ncomx, *gpu_ncomy, *gpu_ncomz;
  // Single molecule arrays
  double *gpu_dx, *gpu_dy, *gpu_dz;
  double **gpu_kx, **gpu_ky, **gpu_kz;
  double **gpu_kxRef, **gpu_kyRef, **gpu_kzRef;
  double **gpu_sumRnew, **gpu_sumInew, **gpu_sumRref, **gpu_sumIref;
  double **gpu_prefact, **gpu_prefactRef;
  double **gpu_hsqr, **gpu_hsqrRef;
  double *gpu_comx, *gpu_comy, *gpu_comz;
  double *gpu_rT11, *gpu_rT12, *gpu_rT13;
  double *gpu_rT22, *gpu_rT23, *gpu_rT33;
  double *gpu_vT11, *gpu_vT12, *gpu_vT13;
  double *gpu_vT22, *gpu_vT23, *gpu_vT33;
  double **gpu_cell_x, **gpu_cell_y, **gpu_cell_z;
  double **gpu_Invcell_x, **gpu_Invcell_y, **gpu_Invcell_z;
  int *gpu_nonOrth;
  double *gpu_aForcex, *gpu_aForcey, *gpu_aForcez;
  double *gpu_mForcex, *gpu_mForcey, *gpu_mForcez;

  double *gpu_t_max;
  double *gpu_r_max;
  double *gpu_BETA;
  double *gpu_mp_coefficient;
  // Currently means translation or rotate MP. Will eventually mean any move.
  int* gpu_move_type;

  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_LJEn;
  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_REn;
  // a flag to prevent translation/rotation of fixed molecules
  int * gpu_moleculeFixed;
  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_coords_x;
  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_coords_y;
  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_coords_z;

  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_com_x;
  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_com_y;
  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_com_z;

  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_aFx;
  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_aFy;
  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_aFz;

  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_mFx;
  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_mFy;
  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_mFz;
  
  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_mTx;
  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_mTy;
  MultiBuffer< DeviceArray<double>, double, buffers > * gpu_mTz;

  double *gpu_aForcex_buffer, *gpu_aForcey_buffer, *gpu_aForcez_buffer;
  double *gpu_mForcex_buffer, *gpu_mForcey_buffer, *gpu_mForcez_buffer;

  double *gpu_mTorquex, *gpu_mTorquey, *gpu_mTorquez;
  int *gpu_inForceRange;

  // Not sure why these need their own arrays
  double *gpu_aForceRecx, *gpu_aForceRecy, *gpu_aForceRecz;
  double *gpu_mForceRecx, *gpu_mForceRecy, *gpu_mForceRecz;

  double *gpu_rMin, *gpu_expConst, *gpu_rMaxSq;

  double *gpu_r_k_x, *gpu_r_k_y, *gpu_r_k_z;
  double *gpu_t_k_x, *gpu_t_k_y, *gpu_t_k_z;

  // Permanent arrays, will try texture after I get global working
  double * gpu_particleCharge;
  int * gpu_particleKind;
  int * gpu_particleMol;

  // lambda structure
  int *gpu_molIndex;
  double *gpu_lambdaVDW, *gpu_lambdaCoulomb;
  bool *gpu_isFraction;

  // To keep molinter working till I port that also
  //int *gpu_cellVector, *gpu_mapParticleToCell, *gpu_cellStartIndex;
  // new pair interaction calculation done on GPU
  MultiBuffer< DeviceArray<int>, int, buffers > * gpu_cellVector;
  MultiBuffer< DeviceArray<int>, int, buffers > * gpu_mapParticleToCell;
  MultiBuffer< DeviceArray<int>, int, buffers > * gpu_cellStartIndex;
  // Fixed as long as volume doesnt change
  // Regenerate after volume moves.
  int *gpu_neighborList;

  // Fixed as long as volume doesnt change
  // Regenerate after volume moves.
  int *gpu_numberOfCells; 
  int *gpu_startOfBoxCellList;
  int *gpu_edgeCells;
  double *gpu_cellSize;
  // Intermediate variable for counting number of molecules in cell
  int *gpu_cellDegrees;
  int *gpu_particleIndices;
  int *gpu_mapParticleToCellSorted;
  int *gpu_Ones;
  int *gpu_CellDegreeSanityCheck;
  int *gpu_IterationsReq;
  size_t zero1 = 0;
  size_t zero2 = 0;

  void    *d_temp_storage_sort_vals;
  void    **d_temp_storage_sort;
  size_t  *temp_storage_bytes_sort;
  
  void    *d_temp_storage_prefix_sum_vals;
  void    **d_temp_storage_prefix_sum;
  size_t  *temp_storage_bytes_prefix_sum;
  // New variables for GPU residence
  double * LJEn, RJEn;
  // For launching kernels
  int cpu_numberOfCells[BOX_TOTAL];
  int cpu_startOfBoxCellList[BOX_TOTAL];

};
#endif
#endif
