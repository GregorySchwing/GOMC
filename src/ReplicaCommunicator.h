#ifndef ReplicaCommunicator_H
#define ReplicaCommunicator_H

#include "GOMC_Config.h"    //For version number

#if GOMC_LIB_MPI
#include <mpi.h>
#endif
#include "XYZArray.h"
#include "Coordinates.h"
#include "COM.h"
#include <cstdint>
#include <cstring>
#include "MoleculeLookup.h"


class ReplicaCommunicator{
  public:
  #if GOMC_LIB_MPI
    ReplicaCommunicator();
    void exchangeXYZArrayNonBlocking(XYZArray * myXYZArray, int exchangePartner);
    void exchangeMolLookupNonBlocking(MoleculeLookup & molLook, int exchangePartner);
    void exchangeMolLookupNonBlocking(uint * molLookup, 
                                      uint molLookupCount, 
                                      uint * boxAndKindStart,
                                      uint boxAndKindStartCount,
                                      uint * boxAndKindSwappable,
                                      uint boxAndKindSwappableCount,
                                      int exchangePartner,
                                      std::vector<uint> & fixedMolecule,
                                      std::vector<uint> & canSwapKind,
                                      std::vector <uint> & canMoveKind,
                                      int * molIndex,
                                      int * molKind,
                                      int * atomIndex,
                                      int * atomKind,
                                      double * atomCharge);
    void exchangeMoleculesNonBlocking(Molecules & mols, int exchangePartner);
    void exchangeMoleculesNonBlocking( uint* start,
                                      uint count,
                                      uint* kIndex,
                                      uint kIndexCount,
                                      uint kindsCount,
                                      uint* countByKind,
                                      char* chain,
                                      double* pairEnCorrections,
                                      double* pairVirCorrections,
                                      int exchangePartner);
  void exchangeBoxDimensionsNonBlocking(BoxDimensions * myBoxDimensions, int exchangePartner);
  #endif
};

#endif