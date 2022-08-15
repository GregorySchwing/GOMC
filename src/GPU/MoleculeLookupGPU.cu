#ifdef GOMC_CUDA
#include "MoleculeLookupGPU.cuh"
#include "cub/cub.cuh"

MoleculeLookupGPU::MoleculeLookupGPU(
                                    uint  & _molLookupCount,
                                    uint  & _atomCount,
                                    uint  & _boxAndKindStartLength,
                                    uint  & _boxAndKindSwappableLength,
                                    uint  & _numKinds,
                                    uint *  _molLookup,  
                                    int  *  _fixedMolecule,                                    
                                    uint *  _boxAndKindStart)
{

    CUMALLOC((void**) gpu_molLookupCount, sizeof(uint));
    CUMALLOC((void**) gpu_atomCount, sizeof(uint));
    CUMALLOC((void**) gpu_boxAndKindStartLength, sizeof(uint));
    CUMALLOC((void**) gpu_boxAndKindSwappableLength, sizeof(uint));
    CUMALLOC((void**) gpu_numKinds, sizeof(uint));

    CUMALLOC((void**) gpu_molLookup, (_molLookupCount + 1) * sizeof(uint));
    CUMALLOC((void**) gpu_fixedMolecule, (_molLookupCount) * sizeof(int));
    CUMALLOC((void**) gpu_boxAndKindStart, (_boxAndKindStartLength) * sizeof(uint));
    CUMALLOC((void**) gpu_mol2Box, _molLookupCount * sizeof(uint));
    CUMALLOC((void**) gpu_numMolsInBox, BOX_TOTAL * sizeof(uint));
    CUMALLOC((void**) gpu_boxMolStartIndex, BOX_TOTAL * sizeof(uint));

    // copy lambda data
    cudaMemcpy(gpu_molLookupCount, &_molLookupCount, sizeof(uint), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_atomCount, &_atomCount, sizeof(uint), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_boxAndKindStartLength, &_boxAndKindStartLength, sizeof(uint), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_boxAndKindSwappableLength, &_boxAndKindSwappableLength, sizeof(uint), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_numKinds, &_numKinds, sizeof(uint), cudaMemcpyHostToDevice);

    cudaMemcpy(gpu_molLookup, _molLookup, (_molLookupCount + 1) * sizeof(uint), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_fixedMolecule, &_fixedMolecule[0], (_molLookupCount) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_boxAndKindStart, _boxAndKindStart, (_boxAndKindStartLength) * sizeof(uint), cudaMemcpyHostToDevice);
    
    uint molBoxStartIndex = 0;
    uint molBoxCount = 0;
    for (int box = 0; box < BOX_TOTAL; ++box){
         molBoxCount = _boxAndKindStart[(box + 1) * _numKinds]
         - _boxAndKindStart[box * _numKinds];
        cudaMemcpy(&gpu_numMolsInBox[box], &molBoxCount, 1 * sizeof(uint), cudaMemcpyHostToDevice);
        cudaMemcpy(&gpu_boxMolStartIndex[box], &molBoxStartIndex, 1 * sizeof(uint), cudaMemcpyHostToDevice);
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