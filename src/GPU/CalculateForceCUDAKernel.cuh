/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.50
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef CALCULATE_FORCE_CUDA_KERNEL
#define CALCULATE_FORCE_CUDA_KERNEL

#ifdef GOMC_CUDA
#include <vector>
#include "XYZArray.h"
#include "BoxDimensions.h"
#include "VariablesCUDA.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#include "CalculateMinImageCUDAKernel.cuh"

using namespace std;

void CallBoxForceGPU(VariablesCUDA *vars,
                     vector<int> cellVector,
                     vector<int> cellStartIndex,
                     std::vector<std::vector<int>> neighborList,
                     XYZArray const &coords,
                     BoxDimensions const &boxAxes,
                     bool electrostatic,
                     vector<double> particleCharge,
                     vector<int> particleKind,
                     vector<int> particleMol,
                     double &REn,
                     double &LJEn,
                     double *aForcex,
                     double *aForcey,
                     double *aForcez,
                     double *mForcex,
                     double *mForcey,
                     double *mForcez,
                     int atomCount,
                     int molCount,
                     bool reset_force,
                     bool copy_back,
                     bool sc_coul,
                     double sc_sigma_6,
                     double sc_alpha,
                     uint sc_power,
                     uint const box);

void CallBoxInterForceGPU(VariablesCUDA *vars,
                          vector<int> cellVector,
                          vector<int> cellStartIndex,
                          std::vector<std::vector<int>> neighborList,
                          XYZArray const &currentCoords,
                          XYZArray const &currentCOM,
                          BoxDimensions const& boxAxes,
                          bool electrostatic,
                          vector<double> &particleCharge,
                          vector<int> &particleKind,
                          vector<int> &particleMol,
                          double &rT11,
                          double &rT12,
                          double &rT13,
                          double &rT22,
                          double &rT23,
                          double &rT33,
                          double &vT11,
                          double &vT12,
                          double &vT13,
                          double &vT22,
                          double &vT23,
                          double &vT33,
                          bool sc_coul,
                          double sc_sigma_6,
                          double sc_alpha,
                          uint sc_power,
                          uint const box);

void CallVirialReciprocalGPU(VariablesCUDA *vars,
                             XYZArray const &currentCoords,
                             XYZArray const &currentCOMDiff,
                             vector<double> &particleCharge,
                             double &rT11,
                             double &rT12,
                             double &rT13,
                             double &rT22,
                             double &rT23,
                             double &rT33,
                             uint imageSize,
                             double constVal,
                             uint box);

__global__ void BoxForceGPU(int *gpu_cellStartIndex,
                            int *gpu_cellVector,
                            int *gpu_neighborList,
                            int numberOfCells,
                            double *gpu_x,
                            double *gpu_y,
                            double *gpu_z,
                            double xAxes,
                            double yAxes,
                            double zAxes,
                            bool electrostatic,
                            double *gpu_particleCharge,
                            int *gpu_particleKind,
                            int *gpu_particleMol,
                            double *gpu_REn,
                            double *gpu_LJEn,
                            int pairSize,
                            double *gpu_sigmaSq,
                            double *gpu_epsilon_Cn,
                            double *gpu_n,
                            int *gpu_VDW_Kind,
                            int *gpu_isMartini,
                            int *gpu_count,
                            double *gpu_rCut,
                            double *gpu_rCutCoulomb,
                            double *gpu_rCutLow,
                            double *gpu_rOn,
                            double *gpu_alpha,
                            int *gpu_ewald,
                            double *gpu_diElectric_1,
                            int *gpu_nonOrth,
                            double *gpu_cell_x,
                            double *gpu_cell_y,
                            double *gpu_cell_z,
                            double *gpu_Invcell_x,
                            double *gpu_Invcell_y,
                            double *gpu_Invcell_z,
                            double *gpu_aForcex,
                            double *gpu_aForcey,
                            double *gpu_aForcez,
                            double *gpu_mForcex,
                            double *gpu_mForcey,
                            double *gpu_mForcez,
                            bool sc_coul,
                            double sc_sigma_6,
                            double sc_alpha,
                            uint sc_power,
                            double *gpu_rMin,
                            double *gpu_rMaxSq,
                            double *gpu_expConst,
                            int box);

__global__ void BoxInterForceGPU(int *gpu_cellStartIndex,
                                 int *gpu_cellVector,
                                 int *gpu_neighborList,
                                 int numberOfCells,
                                 double *gpu_x,
                                 double *gpu_y,
                                 double *gpu_z,
                                 double *gpu_comx,
                                 double *gpu_comy,
                                 double *gpu_comz,
                                 double xAxes,
                                 double yAxes,
                                 double zAxes,
                                 bool electrostatic,
                                 double *gpu_particleCharge,
                                 int *gpu_particleKind,
                                 int *gpu_particleMol,
                                 double *gpu_rT11,
                                 double *gpu_rT12,
                                 double *gpu_rT13,
                                 double *gpu_rT22,
                                 double *gpu_rT23,
                                 double *gpu_rT33,
                                 double *gpu_vT11,
                                 double *gpu_vT12,
                                 double *gpu_vT13,
                                 double *gpu_vT22,
                                 double *gpu_vT23,
                                 double *gpu_vT33,
                                 int pairSize,
                                 double *gpu_sigmaSq,
                                 double *gpu_epsilon_Cn,
                                 double *gpu_n,
                                 int *gpu_VDW_Kind,
                                 int *gpu_isMartini,
                                 int *gpu_count,
                                 double *gpu_rCut,
                                 double *gpu_rCutCoulomb,
                                 double *gpu_rCutLow,
                                 double *gpu_rOn,
                                 double *gpu_alpha,
                                 int *gpu_ewald,
                                 double *gpu_diElectric_1,
                                 double *gpu_cell_x,
                                 double *gpu_cell_y,
                                 double *gpu_cell_z,
                                 double *gpu_Invcell_x,
                                 double *gpu_Invcell_y,
                                 double *gpu_Invcell_z,
                                 int *gpu_nonOrth,
                                 bool sc_coul,
                                 double sc_sigma_6,
                                 double sc_alpha,
                                 uint sc_power,
                                 double *gpu_rMin,
                                 double *gpu_rMaxSq,
                                 double *gpu_expConst,
                                 int box);

__global__ void VirialReciprocalGPU(double *gpu_x,
                                    double *gpu_y,
                                    double *gpu_z,
                                    double *gpu_comDx,
                                    double *gpu_comDy,
                                    double *gpu_comDz,
                                    double *gpu_kxRef,
                                    double *gpu_kyRef,
                                    double *gpu_kzRef,
                                    double *gpu_prefactRef,
                                    double *gpu_hsqrRef,
                                    double *gpu_sumRref,
                                    double *gpu_sumIref,
                                    double *gpu_particleCharge,
                                    double *gpu_rT11,
                                    double *gpu_rT12,
                                    double *gpu_rT13,
                                    double *gpu_rT22,
                                    double *gpu_rT23,
                                    double *gpu_rT33,
                                    double constVal,
                                    uint imageSize,
                                    uint atomNumber);

__device__ double CalcEnForceGPU(double distSq, int kind1, int kind2,
                                 double *gpu_sigmaSq,
                                 double *gpu_n,
                                 double *gpu_epsilon_Cn,
                                 double gpu_rCut,
                                 double gpu_rOn,
                                 int gpu_isMartini,
                                 int gpu_VDW_Kind,
                                 int gpu_count,
                                 double sc_sigma_6,
                                 double sc_alpha,
                                 uint sc_power,
                                 double *gpu_rMin,
                                 double *gpu_rMaxSq,
                                 double *gpu_expConst);

//ElectroStatic Calculation
//**************************************************************//
__device__ double CalcCoulombVirParticleGPU(double distSq, double qi_qj,
    int gpu_ewald, double gpu_alpha,
    int index, double *gpu_sigmaSq,
    bool sc_coul, double sc_sigma_6,
    double sc_alpha, uint sc_power);
__device__ double CalcCoulombVirParticleGPU(double distSq, double qi_qj,
    int gpu_ewald, double gpu_alpha);

//VDW Calculation
//*****************************************************************//
__device__ double CalcVirParticleGPU(double distSq, int index,
                                     double *gpu_sigmaSq, double *gpu_n,
                                     double *gpu_epsilon_Cn,
                                     double sc_sigma_6,
                                     double sc_alpha, uint sc_power);
__device__ double CalcVirParticleGPU(double distSq, int index,
                                     double *gpu_sigmaSq, double *gpu_n,
                                     double *gpu_epsilon_Cn);

// Have to move the implementation for some functions here
// since CUDA doesn't allow __global__ to call __device__
// from different files
// Wanted to call CalcCoulombForceGPU() from CalculateEnergyCUDAKernel.cu file
__device__ inline double CalcCoulombForceGPU(double distSq, double qi_qj,
    int gpu_VDW_Kind, int gpu_ewald,
    int gpu_isMartini,
    double gpu_alpha,
    double gpu_rCutCoulomb,
    double gpu_diElectric_1,
    double *gpu_sigmaSq,
    bool sc_coul,
    double sc_sigma_6,
    double sc_alpha,
    uint sc_power,
    int gpu_count, int kind1,
    int kind2)
{
  if((gpu_rCutCoulomb * gpu_rCutCoulomb) < distSq) {
    return 0.0;
  }

  int index = FlatIndexGPU(kind1, kind2, gpu_count);
  return CalcCoulombVirParticleGPU(distSq, qi_qj, gpu_ewald, gpu_alpha, index,
                                   gpu_sigmaSq, sc_coul, sc_sigma_6, sc_alpha,
                                   sc_power);
}


#endif /*GOMC_CUDA*/
#endif /*CALCULATE_FORCE_CUDA_KERNEL*/
