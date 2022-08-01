#ifndef MOLECULELOOKUP_GPU_H
#define MOLECULELOOKUP_GPU_H
#ifdef GOMC_CUDA

#include "CUDAMemoryManager.cuh"
#include "GOMCEventsProfile.h" // for NVTX profiling
// For cuMemsetD32
//#include <cuda_runtime.h>

class MoleculeLookupGPU {
  public:
    MoleculeLookupGPU(
                                    uint32_t  & _molLookupCount,
                                    uint32_t  & _atomCount,
                                    uint32_t  & _boxAndKindStartLength,
                                    uint32_t  & _boxAndKindSwappableLength,
                                    uint32_t  & _numKinds,
                                    uint32_t *  _molLookup,  
                                    int32_t *  _fixedMolecule,                                    
                                    uint32_t *  _boxAndKindStart);
    ~MoleculeLookupGPU();

  private:
    uint32_t * molLookupCount;
    uint32_t * atomCount;
    uint32_t * boxAndKindStartLength;
    uint32_t * boxAndKindSwappableLength;
    uint32_t * numKinds;

    //array of indices for type Molecule, sorted by box and kind for
    //move selection
    uint32_t* molLookup;

    //index [BOX_TOTAL * kind + box] is the first element of that kind/box in
    //molLookup
    //index [BOX_TOTAL * kind + box + 1] is the element after the end
    //of that kind/box
    uint32_t* boxAndKindStart;
    uint32_t* boxAndKindSwappableCounts;

    int32_t *molIndex; // stores the molecule index for global atom index
    int32_t *atomIndex; // stores the local atom index for global atom index

    int32_t * fixedMolecule; //Molecules that can't move 
    int32_t * canSwapKind; //Kinds that can move intra and inter box
    int32_t * canMoveKind; //Kinds that can move intra box only
};

#endif
#endif