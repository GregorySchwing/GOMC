/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef INTRA_MOL_CELLLIST_H
#define INTRA_MOL_CELLLIST_H
#include "CellList.h"

class BoxDimensions;
class System;

class IntraMolCellList : public CellList {
    public:
    IntraMolCellList(const Molecules& mols,  BoxDimensions& dims, uint molIndex, System& sys);
    bool determineNeedForIMCL();
    System & sysRef;


};

#endif