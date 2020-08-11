#ifdef GOMC_CUDA


#include "PrecisionChecker.cuh"

PrecisionChecker::PrecisionChecker(int i){}

void PrecisionChecker::sortCUDATuples(int * curr, int * neigh, double * val, int numberOfElements){

    thrust::device_vector< int > currVec(curr, curr+numberOfElements);
    thrust::device_vector< int > neighVec(neigh, neigh+numberOfElements);
    thrust::device_vector< double > valVec(val, val+numberOfElements);

    row_vec_dev_cuda_en = currVec;
    col_vec_dev_cuda_en = neighVec;
    val_vec_dev_cuda = valVec;

    A_first_cuda = thrust::make_zip_iterator(thrust::make_tuple(row_vec_dev_cuda_en.begin(), col_vec_dev_cuda_en.begin(), val_vec_dev_cuda.begin()));
    A_last_cuda  = thrust::make_zip_iterator(thrust::make_tuple(row_vec_dev_cuda_en.end(),   col_vec_dev_cuda_en.end(),   val_vec_dev_cuda.end()));

    thrust::sort(A_first_cuda, A_last_cuda, cmpEnergy());

    row_vec_cuda_en = row_vec_dev_cuda_en;
    col_vec_cuda_en = col_vec_dev_cuda_en;
    val_vec_cuda = val_vec_dev_cuda;
}

void PrecisionChecker::sortOMPTuples(int * curr, int * neigh, double * val, int numberOfElements){

    thrust::device_vector< int > currVec(curr, curr+numberOfElements);
    thrust::device_vector< int > neighVec(neigh, neigh+numberOfElements);
    thrust::device_vector< double > valVec(val, val+numberOfElements);

    row_vec_dev_omp_en = currVec;
    col_vec_dev_omp_en = neighVec;
    val_vec_dev_omp = valVec;

    A_first_omp = thrust::make_zip_iterator(thrust::make_tuple(row_vec_dev_omp.begin(), col_vec_dev_omp.begin(), val_vec_dev_omp.begin()));
    A_last_omp  = thrust::make_zip_iterator(thrust::make_tuple(row_vec_dev_omp.end(),   col_vec_dev_omp.end(),   val_vec_dev_omp.end()));

    thrust::sort(A_first_omp, A_last_omp, cmpEnergy());

    row_vec_omp_en = row_vec_dev_omp_en;
    col_vec_omp_en = col_vec_dev_omp_en;
    val_vec_omp = val_vec_dev_omp;
}
/*
void PrecisionChecker::sortCUDATuplesForce(int * curr, int * neigh, double * forceX, double * forceY, double * forceZ, int numberOfElements){

    thrust::device_vector< int > currVec(curr, curr+numberOfElements);
    thrust::device_vector< int > neighVec(neigh, neigh+numberOfElements);
    thrust::device_vector< double > valVecX(forceX, forceX+numberOfElements);
    thrust::device_vector< double > valVecY(forceY, forceY+numberOfElements);
    thrust::device_vector< double > valVecZ(forceZ, forceZ+numberOfElements);

    row_vec_dev_cuda = currVec;
    col_vec_dev_cuda = neighVec;
    valx_vec_dev_cuda = valVecX;
    valy_vec_dev_cuda = valVecY;
    valz_vec_dev_cuda = valVecZ;

    A_first_cuda_force = thrust::make_zip_iterator(thrust::make_tuple(row_vec_dev_cuda.begin(), 
                                                                col_vec_dev_cuda.begin(), 
                                                                valx_vec_dev_cuda.begin(),
                                                                valy_vec_dev_cuda.begin(),
                                                                valz_vec_dev_cuda.begin()));

    A_last_cuda_force  = thrust::make_zip_iterator(thrust::make_tuple(row_vec_dev_cuda.end(),   
                                                                col_vec_dev_cuda.end(),   
                                                                valx_vec_dev_cuda.end(),
                                                                valy_vec_dev_cuda.end(),
                                                                valz_vec_dev_cuda.end()));

    thrust::sort(A_first_cuda_force, A_last_cuda_force, cmpForce());

    row_vec_cuda = row_vec_dev_cuda;
    col_vec_cuda = col_vec_dev_cuda;
    valx_vec_cuda = valx_vec_dev_cuda;
    valy_vec_cuda = valy_vec_dev_cuda;
    valz_vec_cuda = valz_vec_dev_cuda;

}

void PrecisionChecker::sortOMPTuplesForce(int * curr, int * neigh, double * forceX, double * forceY, double * forceZ, int numberOfElements){

    thrust::device_vector< int > currVec(curr, curr+numberOfElements);
    thrust::device_vector< int > neighVec(neigh, neigh+numberOfElements);
    thrust::device_vector< double > valVecX(forceX, forceX+numberOfElements);
    thrust::device_vector< double > valVecY(forceY, forceY+numberOfElements);
    thrust::device_vector< double > valVecZ(forceZ, forceZ+numberOfElements);

    row_vec_dev_omp = currVec;
    col_vec_dev_omp = neighVec;
    valx_vec_dev_omp = valVecX;
    valy_vec_dev_omp = valVecY;
    valz_vec_dev_omp = valVecZ;

    A_first_omp_force = thrust::make_zip_iterator(thrust::make_tuple(row_vec_dev_omp.begin(), 
                                                                col_vec_dev_omp.begin(), 
                                                                valx_vec_dev_omp.begin(),
                                                                valy_vec_dev_omp.begin(),
                                                                valz_vec_dev_omp.begin()));

    A_last_omp_force  = thrust::make_zip_iterator(thrust::make_tuple(row_vec_dev_omp.end(),   
                                                                col_vec_dev_omp.end(),   
                                                                valx_vec_dev_omp.end(),
                                                                valy_vec_dev_omp.end(),
                                                                valz_vec_dev_omp.end()));

    thrust::sort(A_first_omp_force, A_last_omp_force, cmpForce());

    row_vec_omp = row_vec_dev_omp;
    col_vec_omp = col_vec_dev_omp;
    valx_vec_omp = valx_vec_dev_omp;
    valy_vec_omp = valy_vec_dev_omp;
    valz_vec_omp = valz_vec_dev_omp;
}
*/

bool PrecisionChecker::AlmostEqualUlps(float A, float B, int maxUlpsDiff)
{
    Float_t uA(A);
    Float_t uB(B);
    // Different signs means they do not match.
    if (uA.Negative() != uB.Negative())
    {
        // Check for equality to make sure +0==-0
        if (A == B)
            return true;
        return false;
    }
    // Find the difference in ULPs.
    int ulpsDiff = abs(uA.i - uB.i);
    if (ulpsDiff <= maxUlpsDiff)
        return true;
    return false;
}
bool PrecisionChecker::AlmostEqualUlps(double A, double B, int maxUlpsDiff)
{
    Double_t uA(A);
    Double_t uB(B);
    // Different signs means they do not match.
    if (uA.Negative() != uB.Negative())
    {
        // Check for equality to make sure +0==-0
        if (A == B)
            return true;
        return false;
    }
    // Find the difference in ULPs.
    int ulpsDiff = abs(uA.i - uB.i);
    if (ulpsDiff <= maxUlpsDiff)
        return true;
    return false;
}

#endif