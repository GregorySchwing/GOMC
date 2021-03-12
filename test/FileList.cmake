set(TestSources
    test/src/BasicTypesTest.cpp
    test/src/BitLibTest.cpp
    test/src/EndianTest.cpp
    test/src/MolLookupTest.cpp
    test/src/CircuitTester.cpp
    test/src/PSFParserTest.cpp
    test/src/ParallelTemperingTest.cpp
)

set(TestHeaders

)

set(GOMCSources
    src/PDBSetup.cpp
    src/MolSetup.cpp
    src/BondAdjacencyList.cpp
    src/ConfigSetup.cpp
    src/FFSetup.cpp
    src/FFConst.cpp
    src/Reader.cpp
    lib/CircuitFinder.cpp
    src/InputFileReader.cpp
    lib/FloydWarshallCycle.cpp
    src/Simulation.cpp
)

set(GOMCHeaders
    src/PDBSetup.h
    src/MolSetup.h
    src/BondAdjacencyList.h
    src/ConfigSetup.h
    src/FFSetup.h
    src/FFConst.h
    src/Reader.h
    lib/CircuitFinder.h
    src/InputFileReader.h
    lib/FloydWarshallCycle.h
    GOMC_Config.h
    src/Simulation.h
)