#include <gtest/gtest.h>
#include "Simulation.h"
#include<unistd.h> 

#ifdef GOMC_CUDA
TEST(CellListGPU, CheckMETHANOL) {
    int result = chdir("./test/input/Systems/METHANOL_OPLSAA/10K/Standard");
    Simulation base("in_NVT.conf");
    //Simulation base("in_GCMC.conf");
    std::vector<int> cellVector, cellStartIndex, mapParticleToCell;
    std::vector< std::vector<int> > neighborList;
    std::vector<int> cellVectorGPU, cellStartIndexGPU, mapParticleToCellGPU;
    std::vector< std::vector<int> > neighborListGPU;
    uint box = 0;
    base.GetCPUCellList(box,
                        cellVector,
                        cellStartIndex,
                        mapParticleToCell,
                        neighborList);
    std::vector<int> particleIndices;

    base.GetGPUCellList(cellVectorGPU,
                        cellStartIndexGPU,
                        mapParticleToCellGPU,
                        neighborListGPU,
                        particleIndices);
    for (int i = 0; mapParticleToCell.size(); ++i){
        if(mapParticleToCell[i] != mapParticleToCellGPU[i])
            printf("%d %d x\n", mapParticleToCell[i], mapParticleToCellGPU[i]);
        else
            printf("%d %d\n", mapParticleToCell[i], mapParticleToCellGPU[i]);
    
    }
    /*
    for (int i = 0; cellStartIndex.size(); ++i){
        if(cellStartIndex[i] != cellStartIndexGPU[i])
            printf("%d %d\n", cellStartIndex[i], cellStartIndexGPU[i]);
    }
    for (int i = 0; cellVector.size(); ++i){
        if(cellVector[i] != cellVectorGPU[i])
            printf("%d %d\n", cellVector[i], cellVectorGPU[i]);
    }
    */
    EXPECT_EQ(mapParticleToCell, mapParticleToCellGPU);
    EXPECT_EQ(cellStartIndex, cellStartIndexGPU);
    EXPECT_EQ(cellVector, cellVectorGPU);

}
#endif