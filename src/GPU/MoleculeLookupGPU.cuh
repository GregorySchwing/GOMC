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

    uint32_t * gpu_molLookupCount;
    uint32_t * gpu_atomCount;
    uint32_t * gpu_boxMolStartIndex;
    uint32_t * gpu_molsInBox;
    uint32_t * gpu_boxAndKindStartLength;
    uint32_t * gpu_boxAndKindSwappableLength;
    uint32_t * gpu_numKinds;
    uint32_t * gpu_mol2Box;
    //array of indices for type Molecule, sorted by box and kind for
    //move selection
    uint32_t* gpu_molLookup;

    //index [BOX_TOTAL * kind + box] is the first element of that kind/box in
    //molLookup
    //index [BOX_TOTAL * kind + box + 1] is the element after the end
    //of that kind/box
    uint32_t* gpu_boxAndKindStart;
    uint32_t* gpu_boxAndKindSwappableCounts;

    int32_t *gpu_molIndex; // stores the molecule index for global atom index
    int32_t *gpu_atomIndex; // stores the local atom index for global atom index

    int32_t * gpu_fixedMolecule; //Molecules that can't move 
    int32_t * gpu_canSwapKind; //Kinds that can move intra and inter box
    int32_t * gpu_canMoveKind; //Kinds that can move intra box only

    int32_t * gpu_startAtomIdx;

};

#endif
#endif