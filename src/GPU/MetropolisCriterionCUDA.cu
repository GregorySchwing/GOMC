/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at <https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifdef GOMC_CUDA

#include "MetropolisCriterionCUDA.cuh"

void CallBMPAccept(){

}

__global__ void GetCoeffTranslation(   
                            int numberOfMolecules,
                            double t_max,
                            double * BETA,
                            double * t_k,
                            double * molForceRefX,
                            double * molForceRefY,
                            double * molForceRefZ,
                            double * molForceNewX,
                            double * molForceNewY,
                            double * molForceNewZ,
                            double * molForceRecRefX,
                            double * molForceRecRefY,
                            double * molForceRecRefZ,
                            double * molForceRecNewX,
                            double * molForceRecNewY,
                            double * molForceRecNewZ){

    int molNumber = blockIdx.x * blockDim.x + threadIdx.x;
    if (molNumber >= numberOfMolecules) return;
    double t_max4 = t_max*4;
    double w_ratio = 0.0;
    // bf_ = BETA * torque * maxTorque
    double3 bf_old = make_double3   ((molForceRefX[molNumber] + molForceRecRefX[molNumber]),
                                    (molForceRefY[molNumber] + molForceRecRefY[molNumber]),
                                    (molForceRefZ[molNumber] + molForceRecRefZ[molNumber]));
    double3 bf_new = make_double3   ((molForceNewX[molNumber] + molForceRecNewX[molNumber]),
                                    (molForceNewY[molNumber] + molForceRecNewY[molNumber]),
                                    (molForceNewZ[molNumber] + molForceRecNewZ[molNumber]));             
            
    w_ratio += CalculateWRatio(bf_new, bf_old, t_k[molNumber], t_max4) * BETA[0] * t_max;

}

__device__ double CalculateWRatio(  const double3 &lb_new,
                                    const double3 &lb_old,
                                    const double3 &k,
                                    const double max4){
    double w_ratio = 0.0;
    double3 old_var = Subtract(lb_old, k);
    double3 new_var = Add(lb_new, k);

    //Note: we could factor max4 and multiply at the end, but
    //      for the move, where we translate and rotate all molecules,
    //      this method would not work. Hence, I did not factor it.
    // its actually is w_ratio += -1.0* but we simplify it
    w_ratio -= (LengthSq(new_var) / max4);
    // its actually is w_ratio -= -1.0* but we simplify it
    w_ratio += (LengthSq(old_var) / max4);

    return w_ratio;
}


#endif