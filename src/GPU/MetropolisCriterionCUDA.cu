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


void CallGetCoeff(VariablesCUDA *vars,
                    int moveType,
                    int molCount,
                    double * MPCoeff){

    // Zeroes value.
    cudaMemcpy(vars->gpu_mp_coefficient, MPCoeff, 1 * sizeof(double),
              cudaMemcpyHostToDevice);

    int threadsPerBlock = 256;
    int blocksPerGrid = (int)(molCount / threadsPerBlock) + 1;

    if (moveType == MPROTATE){
        BufferAccess<DeviceArray<double>, double, buffers> mTxRef(*(vars->gpu_mTx), curr_state);
        BufferAccess<DeviceArray<double>, double, buffers> mTyRef(*(vars->gpu_mTy), curr_state);
        BufferAccess<DeviceArray<double>, double, buffers> mTzRef(*(vars->gpu_mTz), curr_state);

        BufferAccess<DeviceArray<double>, double, buffers> mTxNew(*(vars->gpu_mTx), next_state);
        BufferAccess<DeviceArray<double>, double, buffers> mTyNew(*(vars->gpu_mTy), next_state);
        BufferAccess<DeviceArray<double>, double, buffers> mTzNew(*(vars->gpu_mTz), next_state);

        GetCoeffRotation<<< blocksPerGrid, threadsPerBlock>>>(molCount,
                                                                vars->gpu_r_max,
                                                                vars->gpu_BETA,
                                                                vars->gpu_mp_coefficient,
                                                                vars->gpu_r_k_x,
                                                                vars->gpu_r_k_y,
                                                                vars->gpu_r_k_z,
                                                                mTxRef->get(),
                                                                mTyRef->get(),
                                                                mTzRef->get(),
                                                                mTxNew->get(),
                                                                mTyNew->get(),
                                                                mTzNew->get());
    } else {
        BufferAccess<DeviceArray<double>, double, buffers> mFxRef(*(vars->gpu_mFx), curr_state);
        BufferAccess<DeviceArray<double>, double, buffers> mFyRef(*(vars->gpu_mFy), curr_state);
        BufferAccess<DeviceArray<double>, double, buffers> mFzRef(*(vars->gpu_mFz), curr_state);

        BufferAccess<DeviceArray<double>, double, buffers> mFxNew(*(vars->gpu_mFx), next_state);
        BufferAccess<DeviceArray<double>, double, buffers> mFyNew(*(vars->gpu_mFy), next_state);
        BufferAccess<DeviceArray<double>, double, buffers> mFzNew(*(vars->gpu_mFz), next_state);

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
    }
    cudaDeviceSynchronize();
    checkLastErrorCUDA(__FILE__, __LINE__);
    cudaMemcpy(MPCoeff, vars->gpu_mp_coefficient, 1 * sizeof(double),
              cudaMemcpyDeviceToHost);
}


__global__ void GetCoeffRotation(   
                            int numberOfMolecules,
                            double * r_max,
                            double * BETA,
                            double * mp_coefficient,
                            double * r_k_x,
                            double * r_k_y,
                            double * r_k_z,
                            double * molTorqueRefX,
                            double * molTorqueRefY,
                            double * molTorqueRefZ,
                            double * molTorqueNewX,
                            double * molTorqueNewY,
                            double * molTorqueNewZ){

    int molNumber = blockIdx.x * blockDim.x + threadIdx.x;
    if (molNumber >= numberOfMolecules) return;

    double r_max4 = r_max[0]*4;
    double w_ratio = 0.0;
    // bf_ = BETA * torque * maxTorque
    double3 bf_old = make_double3   (molTorqueRefX[molNumber] * BETA[0] * r_max[0],
                                    molTorqueRefY[molNumber]  * BETA[0] * r_max[0],
                                    molTorqueRefZ[molNumber]  * BETA[0] * r_max[0]);
    double3 bf_new = make_double3   (molTorqueNewX[molNumber] * BETA[0] * r_max[0],
                                    molTorqueNewY[molNumber]  * BETA[0] * r_max[0],
                                    molTorqueNewZ[molNumber]  * BETA[0] * r_max[0]);             

    double3 k = make_double3   (r_k_x[molNumber],r_k_y[molNumber],r_k_z[molNumber]);

    w_ratio += CalculateWRatio(bf_new, bf_old, k, r_max4);

    atomicAdd(&mp_coefficient[0], w_ratio);

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

    double t_max4 = t_max[0]*4;
    double w_ratio = 0.0;
    // bf_ = BETA * torque * maxTorque
    double3 bf_old = make_double3   ((molForceRefX[molNumber] + molForceRecRefX[molNumber]) * BETA[0] * t_max[0],
                                    (molForceRefY[molNumber] + molForceRecRefY[molNumber])  * BETA[0] * t_max[0],
                                    (molForceRefZ[molNumber] + molForceRecRefZ[molNumber])  * BETA[0] * t_max[0]);
    double3 bf_new = make_double3   ((molForceNewX[molNumber] + molForceRecNewX[molNumber]) * BETA[0] * t_max[0],
                                    (molForceNewY[molNumber] + molForceRecNewY[molNumber])  * BETA[0] * t_max[0],
                                    (molForceNewZ[molNumber] + molForceRecNewZ[molNumber])  * BETA[0] * t_max[0]);             

    double3 k = make_double3   (t_k_x[molNumber],t_k_y[molNumber],t_k_z[molNumber]);

    w_ratio += CalculateWRatio(bf_new, bf_old, k, t_max4);

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