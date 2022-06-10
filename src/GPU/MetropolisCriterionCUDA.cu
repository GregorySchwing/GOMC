/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at <https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifdef GOMC_CUDA

#include "MetropolisCriterionCUDA.cuh"

const int curr_state = 0;
const int next_state = 1;


void CallGetCoeffTranslation(VariablesCUDA *vars,
                    int molCount,
                    double * MPCoeff){

    BufferAccess<DeviceArray<double>, double, buffers> mFxRef(*(vars->gpu_mFx), curr_state);
    BufferAccess<DeviceArray<double>, double, buffers> mFyRef(*(vars->gpu_mFy), curr_state);
    BufferAccess<DeviceArray<double>, double, buffers> mFzRef(*(vars->gpu_mFz), curr_state);

    BufferAccess<DeviceArray<double>, double, buffers> mFxNew(*(vars->gpu_mFx), next_state);
    BufferAccess<DeviceArray<double>, double, buffers> mFyNew(*(vars->gpu_mFy), next_state);
    BufferAccess<DeviceArray<double>, double, buffers> mFzNew(*(vars->gpu_mFz), next_state);

    cudaMemcpy(vars->gpu_mp_coefficient, MPCoeff, 1 * sizeof(double),
              cudaMemcpyHostToDevice);

    int threadsPerBlock = 256;
    int blocksPerGrid = (int)(molCount / threadsPerBlock) + 1;

    GetCoeffTranslation<<< blocksPerGrid, threadsPerBlock>>>(molCount,
                                                            vars->gpu_t_max,
                                                            vars->gpu_BETA,
                                                            vars->gpu_mp_coefficient,
                                                            vars->gpu_t_k_x,
                                                            vars->gpu_t_k_y,
                                                            vars->gpu_t_k_z,
                                                            mFxRef->get(),
                                                            mFyRef->get(),
                                                            mFzRef->get(),
                                                            mFxNew->get(),
                                                            mFyNew->get(),
                                                            mFzNew->get(),
                                                            vars->gpu_mForceRecx,
                                                            vars->gpu_mForceRecy,
                                                            vars->gpu_mForceRecz,
                                                            vars->gpu_mForceRecx,
                                                            vars->gpu_mForceRecy,
                                                            vars->gpu_mForceRecz);
                                                        cudaDeviceSynchronize();
                                                        checkLastErrorCUDA(__FILE__, __LINE__);
    cudaMemcpy(MPCoeff, vars->gpu_mp_coefficient, 1 * sizeof(double),
              cudaMemcpyDeviceToHost);
}

__global__ void GetCoeffTranslation(   
                            int numberOfMolecules,
                            double * t_max,
                            double * BETA,
                            double * mp_coefficient,
                            double * t_k_x,
                            double * t_k_y,
                            double * t_k_z,
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
    printf("Entered method\n");
    double t_max4 = t_max[0]*4;
    double w_ratio = 0.0;
    // bf_ = BETA * torque * maxTorque
    double3 bf_old = make_double3   ((molForceRefX[molNumber] + molForceRecRefX[molNumber]),
                                    (molForceRefY[molNumber] + molForceRecRefY[molNumber]),
                                    (molForceRefZ[molNumber] + molForceRecRefZ[molNumber]));
    double3 bf_new = make_double3   ((molForceNewX[molNumber] + molForceRecNewX[molNumber]),
                                    (molForceNewY[molNumber] + molForceRecNewY[molNumber]),
                                    (molForceNewZ[molNumber] + molForceRecNewZ[molNumber]));             
        printf("Entered bf_new\n");

    double3 k = make_double3   (t_k_x[molNumber],t_k_y[molNumber],t_k_z[molNumber]);
    printf("Entered k\n");

    w_ratio += CalculateWRatio(bf_new, bf_old, k, t_max4) * BETA[0] * t_max[0];
       printf("Entered CalculateWRatio\n");

    atomicAdd(&mp_coefficient[0], w_ratio);

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
//Since the ChangeBuffers and GridAll functions have CPU wrappers, I can't get full GPU residence yet.
/*
__global__ void Accept(   
                            double mp_coefficient,
                            double BETA,
                            double sysPotNew,
                            double sysPotRef){
  double accept = exp(-BETA * (sysPotNew - sysPotRef) + mp_coefficient);
  bool result = true;
  //bool result = prng() < accept;

}
*/
#endif