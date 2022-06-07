/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at <https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef CONSTANT_DEFINITIONS_CUDA_KERNEL
#define CONSTANT_DEFINITIONS_CUDA_KERNEL

#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include "GeomLib.h"
#include "VariablesCUDA.cuh"
#include "EnsemblePreprocessor.h"
#include <cstring>
#include "DoubleBuffer.cuh"
#define GPU_VDW_STD_KIND 0
#define GPU_VDW_SHIFT_KIND 1
#define GPU_VDW_SWITCH_KIND 2
#define GPU_VDW_EXP6_KIND 3
#define MAX_PAIR_SIZE 10000000

void UpdateGPULambda(VariablesCUDA *vars, int *molIndex, double *lambdaVDW,
                    double *lambdaCoulomb, bool *isFraction);
void InitGPUForceField(VariablesCUDA &vars, double const *sigmaSq,
                       double const *epsilon_Cn, double const *n,
                       int VDW_Kind, int isMartini, int count,
                       double Rcut, double const *rCutCoulomb,
                       double RcutLow, double Ron, double const *alpha,
                       int ewald, double diElectric_1);
                       
void InitMPVars(VariablesCUDA *vars,
                std::vector<double> & t_max,
                std::vector<double> & r_max,
                int maxMolNumber);

void InitCoordinatesCUDA(VariablesCUDA *vars, uint atomNumber,
                         double * coords_x,
                         double * coords_y,
                         double * coords_z,
                         double * com_x,
                         double * com_y,
                         double * com_z,
                         uint maxAtomsInMol, uint maxMolNumber,
                         std::vector<double> & particleCharge,
                         std::vector<int> & particleKind,
                         std::vector<int> & particleMol);
void InitEwaldVariablesCUDA(VariablesCUDA *vars, uint imageTotal);
void CopyCurrentToRefCUDA(VariablesCUDA *vars, uint box, uint imageTotal);
void CopyRefToNewCUDA(VariablesCUDA *vars, uint box, uint imageTotal);
void UpdateRecipVecCUDA(VariablesCUDA *vars, uint box);
void UpdateRecipCUDA(VariablesCUDA *vars, uint box);
void UpdateCellBasisCUDA(VariablesCUDA *vars, uint box, double *cellBasis_x,
                         double *cellBasis_y, double *cellBasis_z);
void UpdateInvCellBasisCUDA(VariablesCUDA *vars, uint box,
                            double *invCellBasis_x, double *invCellBasis_y,
                            double *invCellBasis_z);
void DestroyEwaldCUDAVars(VariablesCUDA *vars);
void DestroyCUDAVars(VariablesCUDA *vars);
void InitExp6Variables(VariablesCUDA *vars, double *rMin, double *expConst,
                       double *rMaxSq, uint size);
void InitGPUCellList(VariablesCUDA *vars, 
                    const std::vector<int> &neighborList,
                    const std::vector<int> &numberOfCells,
                    const std::vector<int> &startOfBoxCellList,
                    const std::vector<int> &edgeCells,
                    const std::vector<double> &cellSize);
#endif /*GOMC_CUDA*/
#endif /*CONSTANT_DEFINITIONS_CUDA_KERNEL*/
