#ifndef _replica_state_h
#define _replica_state_h


class SystemPotential;
class Coordinates;
class COM;
class CalcEwald;
class CellList;
class Ewald;
#include "EnergyTypes.h"


class ReplicaState
{
    public:
        ReplicaState(   SystemPotential& potential,
                        Coordinates& coordinates,
                        COM& com,
                        Ewald *calcEwald,
                        CellList& cellList);   
                             
        ReplicaState(   SystemPotential& potential,
                        Coordinates& coordinates,
                        COM& com,
                        Ewald *calcEwald,
                        CellList& cellList):potential(potential), coordinates(coordinates), com(com), calcEwald(calcEwald), cellList(cellList){}
          SystemPotential& potential; //ex
          Coordinates& coordinates; //ex
          COM& com; //ex
          Ewald *calcEwald; //ex
          CellList& cellList; //ex
};

#endif