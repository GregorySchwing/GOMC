#ifdef GOMC_CUDA
#include "MoleculeLookupGPU.cuh"
#include "cub/cub.cuh"

MoleculeLookupGPU::MoleculeLookupGPU(
                                    uint32_t  & _molLookupCount,
                                    uint32_t  & _atomCount,
                                    uint32_t  & _boxAndKindStartLength,
                                    uint32_t  & _boxAndKindSwappableLength,
                                    uint32_t  & _numKinds,
                                    uint32_t *  _molLookup,  
                                    int32_t  *  _fixedMolecule,                                    
                                    uint32_t *  _boxAndKindStart){

    CUMALLOC((void**) molLookupCount, sizeof(uint32_t));
    CUMALLOC((void**) atomCount, sizeof(uint32_t));
    CUMALLOC((void**) boxAndKindStartLength, sizeof(uint32_t));
    CUMALLOC((void**) boxAndKindSwappableLength, sizeof(uint32_t));
    CUMALLOC((void**) numKinds, sizeof(uint32_t));

    CUMALLOC((void**) molLookup, (_molLookupCount + 1) * sizeof(uint32_t));
    CUMALLOC((void**) fixedMolecule, (_molLookupCount) * sizeof(int32_t));
    CUMALLOC((void**) boxAndKindStart, (_boxAndKindStartLength) * sizeof(uint32_t));

    // copy lambda data
    cudaMemcpy(molLookupCount, &_molLookupCount, sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(atomCount, &_atomCount, sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(boxAndKindStartLength, &_boxAndKindStartLength, sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(boxAndKindSwappableLength, &_boxAndKindSwappableLength, sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(numKinds, &_numKinds, sizeof(uint32_t), cudaMemcpyHostToDevice);

    cudaMemcpy(molLookup, _molLookup, (_molLookupCount + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(fixedMolecule, _fixedMolecule, (_molLookupCount) * sizeof(int32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(boxAndKindStart, _boxAndKindStart, (_boxAndKindStartLength) * sizeof(uint32_t), cudaMemcpyHostToDevice);

}

MoleculeLookupGPU::~MoleculeLookupGPU(){
    
    CUFREE(molLookupCount);
    CUFREE(atomCount);
    CUFREE(boxAndKindStartLength);
    CUFREE(boxAndKindSwappableLength);
    CUFREE(numKinds);

    CUFREE(molLookup);
    CUFREE(fixedMolecule);
    CUFREE(boxAndKindStart);

}


#endif