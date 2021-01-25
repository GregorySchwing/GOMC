#include "IntraMolCellList.h"
#include "XYZArray.h"
#include "System.h"

IntraMolCellList::IntraMolCellList(const Molecules& mols,  BoxDimensions& dims, uint molIndex, System& sys):
CellList(mols, dims), sysRef(sys)
{
    GridOne(dims, sysRef.coordinates, sysRef.molLookupRef, 0, molIndex);
    XYZ centerOfCell0 = GetCellZeroCenter(0);
}

