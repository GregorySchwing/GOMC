#include <gtest/gtest.h>
#include "Simulation.h"
#include<unistd.h> 

#ifdef GOMC_CUDA
#if ENSEMBLE == NVT
TEST(CellListGPU, CheckMETHANOL) {
    int result;
    result = chdir("./test/input/Systems/METHANOL_OPLSAA/10K/Standard");
    Simulation base("in_NVT.conf");
    //Simulation base("in_GCMC.conf");
    std::vector<int> cellVector, cellStartIndex, mapParticleToCell;
    std::vector< std::vector<int> > neighborList;
    std::vector<int> cellVectorGPU, cellStartIndexGPU, mapParticleToCellGPU;
    std::vector<int> neighborListGPU;
    uint box = 0;
    base.GetCPUCellList(
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
   /*
    printf("mapParticleToCell.size() %d\n", mapParticleToCell.size());
    printf("mapParticleToCellGPU.size() %d\n", mapParticleToCellGPU.size());
    for (int i = 0; i < mapParticleToCell.size(); ++i){
        if(mapParticleToCell[i] != mapParticleToCellGPU[i])
            printf("mol index %d %d %d x\n", i, mapParticleToCell[i], mapParticleToCellGPU[i]);
        else
            printf("mol index %d %d %d\n", i, mapParticleToCell[i], mapParticleToCellGPU[i]);
    }
 
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
#elif ENSEMBLE == GEMC
TEST(CellListGPU, CheckPEN_HEX) {
    int result;
    result = chdir("./test/input/Systems/PEN_HEX/Base");
    Simulation base("in.conf");
    std::vector<int> cellVector, cellStartIndex, mapParticleToCell;

    std::vector<int> neighborList;

    std::vector<int> cellVectorGPU, cellStartIndexGPU, mapParticleToCellGPU;
    std::vector<int> neighborListGPU;

    base.GetCPUCellList(cellVector,
                        cellStartIndex,
                        mapParticleToCell,
                        neighborList);

    std::vector<int> particleIndices;

    base.GetGPUCellList(cellVectorGPU,
                        cellStartIndexGPU,
                        mapParticleToCellGPU,
                        neighborListGPU,
                        particleIndices);
   /*
    printf("mapParticleToCell.size() %d\n", mapParticleToCell.size());
    printf("mapParticleToCellGPU.size() %d\n", mapParticleToCellGPU.size());
    for (int i = 0; i < mapParticleToCell.size(); ++i){
        if(mapParticleToCell[i] != mapParticleToCellGPU[i])
            printf("mol index %d %d %d x\n", i, mapParticleToCell[i], mapParticleToCellGPU[i]);
        else
            printf("mol index %d %d %d\n", i, mapParticleToCell[i], mapParticleToCellGPU[i]);
    }

    printf("cellStartIndex.size() %d\n", cellStartIndex.size());
    printf("cellStartIndex.size() %d\n", cellStartIndexGPU.size());

    for (int i = 0; i < cellStartIndex.size(); ++i){
        if(cellStartIndex[i] != cellStartIndexGPU[i])
            printf("NOT EQ %d %d\n", cellStartIndex[i], cellStartIndexGPU[i]);
        else 
            printf("EQ %d %d\n", cellStartIndex[i], cellStartIndexGPU[i]);
    }

    for (int i = 0; cellVector.size(); ++i){
        if(cellVector[i] != cellVectorGPU[i])
            printf("%d %d\n", cellVector[i], cellVectorGPU[i]);
    }
    */



    printf("neighborList.size() %d\n", neighborList.size());
    printf("neighborListGPU.size() %d\n", neighborListGPU.size());

    for (int i = 0; i < neighborList.size(); ++i){
        if(neighborList[i] != neighborListGPU[i])
            printf("%d %d\n", neighborList[i], neighborListGPU[i]);
    }
    EXPECT_EQ(mapParticleToCell, mapParticleToCellGPU);
    EXPECT_EQ(cellStartIndex, cellStartIndexGPU);
    EXPECT_EQ(cellVector, cellVectorGPU);
    EXPECT_EQ(neighborList, neighborListGPU);

}
#endif

#endif