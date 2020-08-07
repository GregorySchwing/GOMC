// Author Bruce Dawson
// https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/

#ifndef PRECISION_CHECKER_H
#define PRECISION_CHECKER_H

#include <iostream>
#include <cstdlib>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/functional.h>

using namespace std;

  typedef thrust::tuple<int,int,double> pairTuple;

  typedef thrust::device_vector<int>::iterator                     IntIterator;
  typedef thrust::device_vector<double>::iterator                  DoubleIterator;

  typedef thrust::tuple<IntIterator, IntIterator, DoubleIterator> PairIteratorTuple;
  typedef thrust::zip_iterator<PairIteratorTuple>                   PairIterator;

  struct cmp : public std::binary_function<pairTuple,pairTuple,bool>
  {
      __host__ __device__
          bool operator()(const pairTuple& a, const pairTuple& b) const
          {
              if (thrust::get<0>(a) != thrust::get<0>(b))
                  return thrust::get<0>(a) < thrust::get<0>(b);
              else 
                  return thrust::get<1>(a) < thrust::get<1>(b);
          }
  };

class PrecisionChecker
{
public:
    explicit PrecisionChecker(int i);
    void sortCUDATuples(int * curr, int * neigh, double * val, int numberOfElements);

    union Float_t
    {
        Float_t(float num = 0.0f) : f(num) {}
        // Portable extraction of components.
        bool Negative() const { return i < 0; }
        int32_t RawMantissa() const { return i & ((1 << 23) - 1); }
        int32_t RawExponent() const { return (i >> 23) & 0xFF; }
        int32_t i;
        float f;
    #ifdef _DEBUG
        struct
        {   // Bitfields for exploration. Do not use in production code.
            uint32_t mantissa : 23;
            uint32_t exponent : 8;
            uint32_t sign : 1;
        } parts;
    #endif
    };
    union Double_t
    {
        Double_t(double num = 0.0) : d(num) {}
        // Portable extraction of components.
        bool Negative() const { return i < 0; }
        int64_t RawMantissa() const { return i & (((int64_t)1 << 52) - 1); }
        int64_t RawExponent() const { return (i >> 52) & 0x7FF; }
        // 0x7FF -> 11 1's
        int64_t i;
        double d;
    #ifdef _DEBUG
        struct
        {   // Bitfields for exploration. Do not use in production code.
            uint64_t mantissa : 52;
            uint64_t exponent : 11;
            uint64_t sign : 1;
        } parts;
    #endif
    };

    bool AlmostEqualUlps(float A, float B, int maxUlpsDiff);
    bool AlmostEqualUlps(double A, double B, int maxUlpsDiff);

    thrust::host_vector<int> col_vec_cuda;
    thrust::host_vector<int> row_vec_cuda;
    thrust::host_vector<double> val_vec_cuda;

    thrust::host_vector<int> col_vec_omp;
    thrust::host_vector<int> row_vec_omp;
    thrust::host_vector<double> val_vec_omp;

    thrust::host_vector<int> dimensions_vec;
    thrust::host_vector<int> ones;

    thrust::device_vector<int> col_vec_dev_cuda;
    thrust::device_vector<int> row_vec_dev_cuda;
    thrust::device_vector<double> val_vec_dev_cuda;

    thrust::device_vector<int> col_vec_dev_omp;
    thrust::device_vector<int> row_vec_dev_omp;
    thrust::device_vector<double> val_vec_dev_omp;

// Now we'll create some zip_iterators for A and B
    PairIterator A_first_cuda;   
    PairIterator A_last_cuda;   
    PairIterator A_first_omp;   
    PairIterator A_last_omp;   
};
#endif