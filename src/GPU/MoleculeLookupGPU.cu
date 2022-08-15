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
                                    uint32_t *  _boxAndKindStart)
{

    CUMALLOC((void**) gpu_molLookupCount, sizeof(uint32_t));
    CUMALLOC((void**) gpu_atomCount, sizeof(uint32_t));
    CUMALLOC((void**) gpu_boxAndKindStartLength, sizeof(uint32_t));
    CUMALLOC((void**) gpu_boxAndKindSwappableLength, sizeof(uint32_t));
    CUMALLOC((void**) gpu_numKinds, sizeof(uint32_t));

    CUMALLOC((void**) gpu_molLookup, (_molLookupCount + 1) * sizeof(uint32_t));
    CUMALLOC((void**) gpu_fixedMolecule, (_molLookupCount) * sizeof(int32_t));
    CUMALLOC((void**) gpu_boxAndKindStart, (_boxAndKindStartLength) * sizeof(uint32_t));
    CUMALLOC((void**) gpu_mol2Box, _molLookupCount * sizeof(uint32_t));
    CUMALLOC((void**) gpu_numMolsInBox, BOX_TOTAL * sizeof(uint32_t));
    CUMALLOC((void**) gpu_boxMolStartIndex, BOX_TOTAL * sizeof(uint32_t));

    // copy lambda data
    cudaMemcpy(gpu_molLookupCount, &_molLookupCount, sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_atomCount, &_atomCount, sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_boxAndKindStartLength, &_boxAndKindStartLength, sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_boxAndKindSwappableLength, &_boxAndKindSwappableLength, sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_numKinds, &_numKinds, sizeof(uint32_t), cudaMemcpyHostToDevice);

    cudaMemcpy(gpu_molLookup, _molLookup, (_molLookupCount + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_fixedMolecule, &_fixedMolecule[0], (_molLookupCount) * sizeof(int32_t), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_boxAndKindStart, _boxAndKindStart, (_boxAndKindStartLength) * sizeof(uint32_t), cudaMemcpyHostToDevice);
    
    uint32_t molBoxStartIndex = 0;
    uint32_t molBoxCount = 0;
    for (int box = 0; box < BOX_TOTAL; ++box){
         molBoxCount = _boxAndKindStart[(box + 1) * _numKinds]
         - _boxAndKindStart[box * _numKinds];
        cudaMemcpy(&gpu_numMolsInBox[box], &molBoxCount, 1 * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(&gpu_boxMolStartIndex[box], &molBoxStartIndex, 1 * sizeof(uint32_t), cudaMemcpyHostToDevice);
        molBoxStartIndex += molBoxCount;
    }

}

MoleculeLookupGPU::~MoleculeLookupGPU(){
    
    CUFREE(gpu_molLookupCount);
    CUFREE(gpu_atomCount);
    CUFREE(gpu_boxAndKindStartLength);
    CUFREE(gpu_boxAndKindSwappableLength);
    CUFREE(gpu_numKinds);
    CUFREE(gpu_numMolsInBox);
    CUFREE(gpu_boxMolStartIndex);

    CUFREE(gpu_mol2Box);
    CUFREE(gpu_molLookup);
    CUFREE(gpu_fixedMolecule);
    CUFREE(gpu_boxAndKindStart);

}


#endif