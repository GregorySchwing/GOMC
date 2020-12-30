#include <gtest/gtest.h>

#define MOLCOUNT 20
#define KINDCOUNT 2
#define BOXCOUNT 2
#define RATIO_SWAPPABLE 2

/* So we can generate the underlying arrays representing the molecules */
#define private public
#include "MoleculeLookup.h"
#undef private

TEST(MolLookupTest, CheckShiftMolTest) {

    MoleculeLookup * ml = new MoleculeLookup();

    ml->numKinds = KINDCOUNT;

    ml->fixedAtom = std::vector<uint>(MOLCOUNT);
    ml->boxAndKindStart = new uint[ml->numKinds * BOXCOUNT + 1];
    ml->molLookup = new uint[MOLCOUNT];
    ml->molLookupCount = MOLCOUNT;
    ml->boxAndKindSwappableCounts = new uint[ml->numKinds * BOXCOUNT];

    for (int i = 0; i < MOLCOUNT; i++){
        ml->molLookup[i] = i;
        if (i % (MOLCOUNT / BOXCOUNT / KINDCOUNT) < 
                MOLCOUNT / BOXCOUNT / KINDCOUNT / RATIO_SWAPPABLE)
            ml->fixedAtom[i] = 0;
        else
            ml->fixedAtom[i] = 2;
        //std::cout << ml->fixedAtom[i] << std::endl;
    }

    std::vector<uint> copyOfFixedAtom = ml->fixedAtom;


    for (int i = 0; i < ml->numKinds * BOXCOUNT + 1; i++){
        ml->boxAndKindStart[i] = i * (MOLCOUNT / BOXCOUNT / KINDCOUNT);
        //std::cout << ml->boxAndKindStart[i] << std::endl;
    }

    for (int i = 0; i < ml->numKinds * BOXCOUNT; i++){
        ml->boxAndKindSwappableCounts[i] = 
            (MOLCOUNT / BOXCOUNT / KINDCOUNT / RATIO_SWAPPABLE);
        //std::cout << ml->boxAndKindSwappableCounts[i] << std::endl;
    }

    uint mol =  0;
    uint kind = 0;
    uint currentBox = 0;
    uint intoBox = 1;

    uint pos = 0;

    uint temp = copyOfFixedAtom[mol];
    uint startIndex = kind + intoBox * ml->numKinds;

    copyOfFixedAtom.insert(copyOfFixedAtom.begin() + 
        ml->boxAndKindStart[startIndex] + 
        ml->boxAndKindSwappableCounts[startIndex], temp);
    
    if (mol < ml->boxAndKindStart[startIndex] + 
        ml->boxAndKindSwappableCounts[startIndex]){
        copyOfFixedAtom.erase(copyOfFixedAtom.begin() + pos);
    } else {
        copyOfFixedAtom.erase(copyOfFixedAtom.begin() + pos + 1);
    }

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->molLookup[i] << " ";
    }
    std::cout << std::endl;
    
    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->fixedAtom[ml->molLookup[i]] << " ";
    }
    std::cout << std::endl;

    std::cout << "Shifting " << mol << std::endl;
    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->molLookup[i] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->fixedAtom[ml->molLookup[i]] << " ";
    }
    std::cout << std::endl;

    ml->ShiftMolBox(mol, currentBox, intoBox, kind);

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->molLookup[i] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->fixedAtom[ml->molLookup[i]] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        EXPECT_EQ(ml->fixedAtom[ml->molLookup[i]], copyOfFixedAtom[i]);
    }
        copyOfFixedAtom.insert(copyOfFixedAtom.begin() + 
        ml->boxAndKindStart[startIndex] + 
        ml->boxAndKindSwappableCounts[startIndex], temp);
    
    if (mol < ml->boxAndKindStart[startIndex] + 
        ml->boxAndKindSwappableCounts[startIndex]){
        copyOfFixedAtom.erase(copyOfFixedAtom.begin() + pos);
    } else {
        copyOfFixedAtom.erase(copyOfFixedAtom.begin() + pos + 1);
    }

    mol =  1;
    std::cout << "Shifting " << mol << std::endl;
    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->molLookup[i] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->fixedAtom[ml->molLookup[i]] << " ";
    }
    std::cout << std::endl;

    ml->ShiftMolBox(mol, currentBox, intoBox, kind);

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->molLookup[i] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->fixedAtom[ml->molLookup[i]] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        EXPECT_EQ(ml->fixedAtom[ml->molLookup[i]], copyOfFixedAtom[i]);
    }


    mol =  16;
    kind = 1;
    currentBox = 1;
    intoBox = 0;

    pos = 16;

    startIndex = kind + intoBox * ml->numKinds;

    temp = 0;

    for (int i = 0; i < MOLCOUNT; i++){
        EXPECT_EQ(ml->fixedAtom[ml->molLookup[i]], copyOfFixedAtom[i]);
    }

    copyOfFixedAtom.insert(copyOfFixedAtom.begin() + 
    ml->boxAndKindStart[startIndex] + 
    ml->boxAndKindSwappableCounts[startIndex], temp);
    
    if (pos < ml->boxAndKindStart[startIndex] + 
        ml->boxAndKindSwappableCounts[startIndex]){
        copyOfFixedAtom.erase(copyOfFixedAtom.begin() + pos);
    } else {
        copyOfFixedAtom.erase(copyOfFixedAtom.begin() + pos + 1);
    }
    std::cout << "Shifting " << mol << std::endl;
    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->molLookup[i] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->fixedAtom[ml->molLookup[i]] << " ";
    }
    std::cout << std::endl;

    ml->ShiftMolBox(mol, currentBox, intoBox, kind);

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->molLookup[i] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->fixedAtom[ml->molLookup[i]] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        EXPECT_EQ(ml->fixedAtom[ml->molLookup[i]], copyOfFixedAtom[i]);
    }

    pos =  16;

    copyOfFixedAtom.insert(copyOfFixedAtom.begin() + 
    ml->boxAndKindStart[startIndex] + 
    ml->boxAndKindSwappableCounts[startIndex], temp);
    
    if (pos < ml->boxAndKindStart[startIndex] + 
        ml->boxAndKindSwappableCounts[startIndex]){
        copyOfFixedAtom.erase(copyOfFixedAtom.begin() + pos);
    } else {
        copyOfFixedAtom.erase(copyOfFixedAtom.begin() + pos + 1);
    }

    mol =  15;
    std::cout << "Shifting " << mol << std::endl;
    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->molLookup[i] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->fixedAtom[ml->molLookup[i]] << " ";
    }
    std::cout << std::endl;

    ml->ShiftMolBox(mol, currentBox, intoBox, kind);

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->molLookup[i] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->fixedAtom[ml->molLookup[i]] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        EXPECT_EQ(ml->fixedAtom[ml->molLookup[i]], copyOfFixedAtom[i]);
    }

    mol =  11;
    kind = 0;
    currentBox = 1;
    intoBox = 0;
    startIndex = kind + intoBox * ml->numKinds;

    copyOfFixedAtom.insert(copyOfFixedAtom.begin() + 
    ml->boxAndKindStart[startIndex] + 
    ml->boxAndKindSwappableCounts[startIndex], temp);
    
    if (mol < ml->boxAndKindStart[startIndex] + 
        ml->boxAndKindSwappableCounts[startIndex]){
        copyOfFixedAtom.erase(copyOfFixedAtom.begin() + mol);
    } else {
        copyOfFixedAtom.erase(copyOfFixedAtom.begin() + mol + 1);
    }

    mol =  1;
    std::cout << "Shifting " << mol << std::endl;
    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->molLookup[i] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->fixedAtom[ml->molLookup[i]] << " ";
    }
    std::cout << std::endl;

    ml->ShiftMolBox(mol, currentBox, intoBox, kind);

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->molLookup[i] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        if (i == ml->NumInBox(0))
            std::cout << " - ";
        std::cout << ml->fixedAtom[ml->molLookup[i]] << " ";
    }
    std::cout << std::endl;

    for (int i = 0; i < MOLCOUNT; i++){
        EXPECT_EQ(ml->fixedAtom[ml->molLookup[i]], copyOfFixedAtom[i]);
    }
}
