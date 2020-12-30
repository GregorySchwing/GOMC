set(TestSources
    test/src/BasicTypesTest.cpp
    test/src/BitLibTest.cpp
    test/src/PSFParserTest.cpp
    test/src/MolLookupTest.cpp
)

set(TestGOMCSources
    src/PDBSetup.cpp
    src/MolSetup.cpp
    src/BondAdjacencyList.cpp
    src/FFSetup.cpp
    src/FFConst.cpp
    src/Reader.cpp
    src/MoleculeLookup.cpp
)

set(TestHeaders
    src/PDBSetup.h
    src/MolSetup.h
    src/BondAdjacencyList.h
    src/ConfigSetup.h
    src/FFSetup.h
    src/FFConst.h
    src/Reader.h
    src/MoleculeLookup.h
)