#ifndef MOLECULELOOKUP_GPU_H
#define MOLECULELOOKUP_GPU_H
#ifdef GOMC_CUDA

#include "CUDAMemoryManager.cuh"
#include "EnsemblePreprocessor.h"
#include "GOMCEventsProfile.h" // for NVTX profiling
// For cuMemsetD32
//#include <cuda_runtime.h>

class MoleculeLookupGPU {
  public:
    MoleculeLookupGPU(
                                    uint  & _molLookupCount,
                                    uint  & _atomCount,
                                    uint  & _boxAndKindStartLength,
                                    uint  & _boxAndKindSwappableLength,
                                    uint  & _numKinds,
                                    uint *  _molLookup,  
                                    int *  _fixedMolecule,                                    
                                    uint *  _boxAndKindStart);
    ~MoleculeLookupGPU();

  private:

    uint * gpu_molLookupCount;
    uint * gpu_atomCount;
    uint * gpu_boxMolStartIndex;
    uint * gpu_numMolsInBox;
    uint * gpu_boxAndKindStartLength;
    uint * gpu_boxAndKindSwappableLength;
    uint * gpu_numKinds;
    uint * gpu_mol2Box;
    //array of indices for type Molecule, sorted by box and kind for
    //move selection
    uint* gpu_molLookup;

    //index [BOX_TOTAL * kind + box] is the first element of that kind/box in
    //molLookup
    //index [BOX_TOTAL * kind + box + 1] is the element after the end
    //of that kind/box
    uint* gpu_boxAndKindStart;
    uint* gpu_boxAndKindSwappableCounts;

    int *gpu_molIndex; // stores the molecule index for global atom index
    int *gpu_atomIndex; // stores the local atom index for global atom index

    int * gpu_fixedMolecule; //Molecules that can't move 
    int * gpu_canSwapKind; //Kinds that can move intra and inter box
    int * gpu_canMoveKind; //Kinds that can move intra box only

    int * gpu_startAtomIdx;
    friend class CellListGPU;
};

#endif
#endif