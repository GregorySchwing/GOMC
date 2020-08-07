#include "PrecisionChecker.cuh"

PrecisionChecker::PrecisionChecker(int i){}

void PrecisionChecker::sortCUDATuples(int * curr, int * neigh, double * val, int numberOfElements){

    thrust::device_vector< int > currVec(curr, curr+numberOfElements);
    thrust::device_vector< int > neighVec(neigh, neigh+numberOfElements);
    thrust::device_vector< double > valVec(val, val+numberOfElements);

    row_vec_dev_cuda = currVec;
    col_vec_dev_cuda = neighVec;
    val_vec_dev_cuda = valVec;

    A_first_cuda = thrust::make_zip_iterator(thrust::make_tuple(row_vec_dev_cuda.begin(), col_vec_dev_cuda.begin(), val_vec_dev_cuda.begin()));
    A_last_cuda  = thrust::make_zip_iterator(thrust::make_tuple(row_vec_dev_cuda.end(),   col_vec_dev_cuda.end(),   val_vec_dev_cuda.end()));

    thrust::sort(A_first_cuda, A_last_cuda, cmp());

    row_vec_cuda = row_vec_dev_cuda;
    col_vec_cuda = col_vec_dev_cuda;
    val_vec_cuda = val_vec_dev_cuda;
}

void PrecisionChecker::sortOMPTuples(int * curr, int * neigh, double * val, int numberOfElements){

    thrust::device_vector< int > currVec(curr, curr+numberOfElements);
    thrust::device_vector< int > neighVec(neigh, neigh+numberOfElements);
    thrust::device_vector< double > valVec(val, val+numberOfElements);

    row_vec_dev_omp = currVec;
    col_vec_dev_omp = neighVec;
    val_vec_dev_omp = valVec;

    A_first_omp = thrust::make_zip_iterator(thrust::make_tuple(row_vec_dev_omp.begin(), col_vec_dev_omp.begin(), val_vec_dev_omp.begin()));
    A_last_omp  = thrust::make_zip_iterator(thrust::make_tuple(row_vec_dev_omp.end(),   col_vec_dev_omp.end(),   val_vec_dev_omp.end()));

    thrust::sort(A_first_omp, A_last_omp, cmp());

    row_vec_omp = row_vec_dev_omp;
    col_vec_omp = col_vec_dev_omp;
    val_vec_omp = val_vec_dev_omp;
}


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