#include "ReplicaCommunicator.h"

#if GOMC_LIB_MPI
ReplicaCommunicator::ReplicaCommunicator(){}

void ReplicaCommunicator::exchangeXYZArrayNonBlocking(XYZArray * myXYZArray, int exchangePartner)
{
    XYZArray outBuffer(*myXYZArray);
    MPI_Request mpi_req;

    int myXYZArrayCount = outBuffer.Count();
    int otherXYZCount;

    MPI_Isend(&myXYZArrayCount, 1 * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(&otherXYZCount, 1 * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    XYZArray inBuffer(otherXYZCount);

    MPI_Isend(outBuffer.x, outBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBuffer.x, inBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(outBuffer.y, outBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBuffer.y, inBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(outBuffer.z, outBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBuffer.z, inBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    swap(*myXYZArray, inBuffer);
}

void ReplicaCommunicator::exchangeMolLookupNonBlocking(MoleculeLookup & molLook, int exchangePartner){
    exchangeMolLookupNonBlocking(molLook.molLookup, 
                                molLook.molLookupCount, 
                                molLook.boxAndKindStart,
                                molLook.boxAndKindStartCount,
                                molLook.boxAndKindSwappable,
                                molLook.boxAndKindSwappableCount,
                                exchangePartner,
                                molLook.fixedMolecule,
                                molLook.canSwapKind,
                                molLook.canMoveKind,
                                molLook.molIndex,
                                molLook.molKind,
                                molLook.atomIndex,
                                molLook.atomKind,
                                molLook.atomCharge);
}

/* Once I convert checkpoint from modifying output to reproducing the original
   Molecules and Molecule Lookup, I will only need to swap the molLookup, boxAndKindStart,
   and boxAndKindSwappable
*/

void ReplicaCommunicator::exchangeMolLookupNonBlocking(uint * molLookup, 
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
                                                        double * atomCharge)

{
    int atomCount = fixedMolecule.size();

    uint * outBufferMolLookup = new uint[molLookupCount];
    uint * outBufferBoxAndKindStart = new uint[boxAndKindStartCount];
    uint * outBufferBoxAndKindSwappable = new uint[boxAndKindSwappableCount];
    std::vector<uint> outBufferFixedMolecule(fixedMolecule);
    std::vector<uint> outBufferCanSwapKind(canSwapKind);
    std::vector<uint> outBufferCanMoveKind(canMoveKind);
    int * outBufferMolIndex = new int[molLookupCount];
    int * outBufferMolKind = new int[molLookupCount];
    int * outBufferAtomIndex = new int[atomCount];
    int * outBufferAtomKind = new int[atomCount];
    int * outBufferAtomCharge = new int[atomCount];

    std::memcpy(outBufferMolLookup, molLookup, molLookupCount);
    std::memcpy(outBufferBoxAndKindStart, boxAndKindStart, boxAndKindStartCount);
    std::memcpy(outBufferBoxAndKindSwappable, boxAndKindSwappable, boxAndKindSwappableCount);
    std::memcpy(outBufferMolIndex, molIndex, molLookupCount);
    std::memcpy(outBufferMolKind, molKind, molLookupCount);
    std::memcpy(outBufferAtomIndex, atomIndex, atomCount);
    std::memcpy(outBufferAtomKind, atomKind, atomCount);
    std::memcpy(outBufferAtomCharge, atomCharge, atomCount);

    uint * inBufferMolLookup = new uint[molLookupCount];
    uint * inBufferBoxAndKindStart = new uint[boxAndKindStartCount];
    uint * inBufferBoxAndKindSwappable = new uint[boxAndKindSwappableCount];
    std::vector<uint> inBufferFixedMolecule(fixedMolecule.size());
    std::vector<uint> inBufferCanSwapKind(canSwapKind.size());
    std::vector<uint> inBufferCanMoveKind(canMoveKind.size());
    int * inBufferMolIndex = new int[molLookupCount];
    int * inBufferMolKind = new int[molLookupCount];
    int * inBufferAtomIndex = new int[atomCount];
    int * inBufferAtomKind = new int[atomCount];
    double * inBufferAtomCharge = new double[atomCount];

    MPI_Request mpi_req;

    MPI_Isend(outBufferMolLookup, molLookupCount * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBufferMolLookup, molLookupCount * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(outBufferBoxAndKindStart, boxAndKindStartCount * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBufferBoxAndKindStart, boxAndKindStartCount * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(outBufferBoxAndKindSwappable, boxAndKindSwappableCount * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBufferBoxAndKindStart, boxAndKindSwappableCount * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(&outBufferFixedMolecule.front(), outBufferFixedMolecule.size() * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(&inBufferFixedMolecule.front(), inBufferFixedMolecule.size() * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(&outBufferCanSwapKind.front(), outBufferCanSwapKind.size() * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(&inBufferCanSwapKind.front(), inBufferCanSwapKind.size() * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(&outBufferCanMoveKind.front(), outBufferCanMoveKind.size() * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(&inBufferCanMoveKind.front(), inBufferCanMoveKind.size() * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
    
    MPI_Isend(outBufferMolIndex, molLookupCount * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBufferMolIndex, molLookupCount * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(outBufferMolKind, molLookupCount * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBufferMolKind, molLookupCount * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(outBufferAtomIndex, atomCount * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBufferAtomIndex, atomCount * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(outBufferAtomKind, atomCount * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBufferAtomKind, atomCount * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(outBufferAtomCharge, atomCount * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBufferAtomCharge, atomCount * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    std::memcpy(molLookup, inBufferMolLookup, molLookupCount);
    std::memcpy(boxAndKindStart, inBufferBoxAndKindStart, boxAndKindStartCount);
    std::memcpy(boxAndKindSwappable, inBufferBoxAndKindSwappable, boxAndKindStartCount);

    fixedMolecule = inBufferFixedMolecule;
    canSwapKind = inBufferCanSwapKind;
    canMoveKind = inBufferCanMoveKind;

    std::memcpy(molIndex, inBufferMolIndex, boxAndKindStartCount);
    std::memcpy(molKind, inBufferMolKind, boxAndKindStartCount);
    std::memcpy(atomIndex, inBufferAtomIndex, atomCount);
    std::memcpy(atomKind, inBufferAtomKind, atomCount);
    std::memcpy(atomCharge, inBufferAtomCharge, atomCount);

    delete outBufferMolLookup;
    delete outBufferBoxAndKindStart;
    delete outBufferBoxAndKindSwappable;
    delete outBufferMolIndex;
    delete outBufferMolKind;
    delete outBufferAtomIndex;
    delete outBufferAtomKind;
    delete outBufferAtomCharge;

    delete inBufferMolLookup;
    delete inBufferBoxAndKindStart;
    delete inBufferBoxAndKindSwappable;
    delete inBufferMolIndex;
    delete inBufferMolKind;
    delete inBufferAtomIndex;
    delete inBufferAtomKind;
    delete inBufferAtomCharge;

}

/* Potential problems are redefinition of kind 0, kind 1, dependending on replica inputs */
void ReplicaCommunicator::exchangeMoleculesNonBlocking(Molecules & mols, int exchangePartner){

}
void ReplicaCommunicator::exchangeMoleculesNonBlocking( uint* start,
                                                        uint count,
                                                        uint* kIndex,
                                                        uint kIndexCount,
                                                        uint kindsCount,
                                                        uint* countByKind,
                                                        char* chain,
                                                        double* pairEnCorrections,
                                                        double* pairVirCorrections,
                                                        int exchangePartner){
    int atomCount = start[count];

    uint * outBufferStart = new uint[count];
    uint * outBufferKIndex = new uint[kIndexCount];
    uint * outBufferCountByKind = new uint[kindsCount];
    char * outBufferChain = new char[atomCount];
    double * outBufferPairEnCorrections[kindsCount * kindsCount];
    double * outBufferPairVirCorrections[kindsCount * kindsCount];

}

void ReplicaCommunicator::exchangeBoxDimensionsNonBlocking(BoxDimensions * myXYZArray, int exchangePartner)
{
    /*
    XYZArray outBuffer(*myXYZArray);
    MPI_Request mpi_req;

    int myXYZArrayCount = outBuffer.Count();
    int otherXYZCount;

    MPI_Isend(&myXYZArrayCount, 1 * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(&otherXYZCount, 1 * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    XYZArray inBuffer(otherXYZCount);

    MPI_Isend(outBuffer.x, outBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBuffer.x, inBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(outBuffer.y, outBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBuffer.y, inBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(outBuffer.z, outBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBuffer.z, inBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    swap(*myXYZArray, inBuffer);
    */
}

#endif