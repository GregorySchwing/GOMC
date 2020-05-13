#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<thrust/device_ptr.h>
#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include<thrust/sort.h>
#include <ctime>
#include<curand_kernel.h>
#include<curand.h>
#include<random>
using namespace std;
#define numOfArrays 10000
#define maxElements 1000
int main (){
    const int range_from = 0;
    const unsigned int range_to = 2147483647; //2^31 - 1
    random_device rand_dev;
    mt19937 generator(rand_dev());
    uniform_int_distribution<int> distr(range_from, range_to);
    thrust::host_vector<int> h_vec(numOfArrays*maxElements);
    thrust::device_vector<int> d_vec;
    thrust::host_vector<int> h_keys(numOfArrays*maxElements);
    thrust::device_vector<int> d_keys;
    srand(time(NULL));
    size_t f, t;
    cudaSetDevice(0);
    cudaMemGetInfo(&f, &t);
    //new data gens
    for(int i = 0; i < numOfArrays; i++){
        for(int j = 0; j < maxElements; j++){
            h_vec[i*maxElements+j] = distr(generator) ;
        }
    }
                
        //initializing the keys
    int timeKeys = clock();
    for(int i = 0; i < numOfArrays; i++){
    for(int j = 0; j < maxElements; j++){
        h_keys[i*maxElements+j] = i;
        }
    }
    timeKeys = clock()-timeKeys;
    //copying the data to device
    d_vec = h_vec;
    d_keys = h_keys;
    int start_s=clock();
    thrust::stable_sort(d_keys.begin(), d_keys.end());
    int stop_s=clock(); 
    //copying back
    h_vec = d_vec;
    cout  << ((stop_s-start_s)+timeKeys)/double(CLOCKS_PER_SEC)*1000 << endl;
    unsigned* my_device_pointer = thrust::raw_pointer_cast(&d_keys[0]);
    return 0;
}