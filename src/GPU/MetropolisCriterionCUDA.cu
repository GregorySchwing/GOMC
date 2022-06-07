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
                            double * t_max,
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
                            double * molForceRecNewZ
                        ){
          //forceReal =  make_double3(0.0, 0.0, 0.0);

}

__device__ double CalculateWRatio(  const double3 &lb_new,
                                    const double3 &lb_old,
                                    const double3 &k,
                                    const double max4){
    double w_ratio = 0.0;
    //double3 old_var = lb_old - k;
    //double3 new_var = lb_new + k;

    //Note: we could factor max4 and multiply at the end, but
    //      for the move, where we translate and rotate all molecules,
    //      this method would not work. Hence, I did not factor it.
    // its actually is w_ratio += -1.0* but we simplify it
    //w_ratio -= (LengthSq(new_var) / max4);
    // its actually is w_ratio -= -1.0* but we simplify it
    //w_ratio += (LengthSq(old_var) / max4);

    return w_ratio;
}


#endif