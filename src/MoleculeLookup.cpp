/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "MoleculeLookup.h" //Header spec.
#include "EnsemblePreprocessor.h" //For box total
#include "PRNG.h" //For index selection
#include "PDBSetup.h" //For init.
#include "Molecules.h" //For init.
#include <algorithm>
#include <utility>
#include <iostream>

void MoleculeLookup::Init(const Molecules& mols,
                          const pdb_setup::Atoms& atomData)
{
  numKinds = mols.GetKindsCount();
  molLookup = new uint[mols.count];
  molLookupCount = mols.count;

  //+1 to store end value
  boxAndKindStart = new uint[numKinds * BOX_TOTAL + 1];
  boxAndKindStartCount = numKinds * BOX_TOTAL + 1;

  /* Assume the swappable molecules are sorted 
  at the beggining of each boxAndKindStart section.

  AR Box 1
  a b c d e f g
  ^ ^ ^
  Swappable
  boxAndKindSwappableCounts[BOX_TOTAL * kind + box] is 
  the number of swappable elements of that kind/box in molLookup

  boxAndKindSwappableCounts[0] = 3
  
  If we keep the swappable ones sorted at the front,
  then all we need is counts for the swappable ones. */
  boxAndKindSwappableCounts = new uint[numKinds * BOX_TOTAL];
  for (int i = 0; i < numKinds * BOX_TOTAL; i++){
    boxAndKindSwappableCounts[i] = 0;
  }

  // vector[box][kind] = list of mol indices for kind in box
  std::vector<std::vector<std::vector<uint> > > indexVector;
  indexVector.resize(BOX_TOTAL);
  fixedAtom.resize(mols.count);

  for (uint b = 0; b < BOX_TOTAL; ++b) {
    indexVector[b].resize(numKinds);
  }

  for(uint m = 0; m < mols.count; ++m) {
    uint box = atomData.box[atomData.startIdxRes[m]];
    uint kind = mols.kIndex[m];

    fixedAtom[m] = atomData.molBeta[m];
    /* increment count if swappable */
    if (fixedAtom[m] == 0){
      //indexVector[box][kind].insert(indexVector[box][kind].begin(), m);
      indexVector[box][kind].push_back(m);
      boxAndKindSwappableCounts[box * numKinds + kind]++;
    } else {
      indexVector[box][kind].push_back(m);
    }

    //Find the kind that can be swap(beta == 0) or move(beta == 0 or 2)
    if(fixedAtom[m] == 0) {
      if(std::find(canSwapKind.begin(), canSwapKind.end(), kind) ==
          canSwapKind.end())
        canSwapKind.push_back(kind);

      if(std::find(canMoveKind.begin(), canMoveKind.end(), kind) ==
          canMoveKind.end())
        canMoveKind.push_back(kind);

    } else if(fixedAtom[m] == 2) {
      if(std::find(canMoveKind.begin(), canMoveKind.end(), kind) ==
          canMoveKind.end())
        canMoveKind.push_back(kind);

    }
  }

  for (uint b = 0; b < BOX_TOTAL; ++b) {
    for (uint k = 0; k < numKinds; ++k) {
      std::vector<std::pair<uint, uint>> pairsForSorting;
      for(int i = 0; i < indexVector[b][k].size(); i++){
        pairsForSorting.push_back(std::make_pair(fixedAtom[indexVector[b][k][i]], indexVector[b][k][i]));
      }
      std::sort(pairsForSorting.begin(), pairsForSorting.end());
      for(int i = 0; i < indexVector[b][k].size(); i++){
        indexVector[b][k][i] = pairsForSorting[i].second;
      }
    }
  }


/* We need to initialize a boxAndKindStart for the swappable counts */

  uint* progress = molLookup;
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    for (uint k = 0; k < numKinds; ++k) {
      boxAndKindStart[b * numKinds + k] = progress - molLookup;
      progress = std::copy(indexVector[b][k].begin(),
                           indexVector[b][k].end(), progress);
    }
  }
  boxAndKindStart[numKinds * BOX_TOTAL] = mols.count;
}

uint MoleculeLookup::NumInBox(const uint box) const
{
  return boxAndKindStart[(box + 1) * numKinds]
         - boxAndKindStart[box * numKinds];
}

void MoleculeLookup::TotalAndDensity
(uint * numByBox, uint * numByKindBox, double * molFractionByKindBox,
 double * densityByKindBox, double const*const volInv) const
{
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    numByBox[b] = 0.0;
    for (uint k = 0; k < numKinds; ++k) {
      uint numMK = NumKindInBox(k, b);
      uint mkIdx = k + numKinds * b;
      numByKindBox[mkIdx] = numMK;
      densityByKindBox[mkIdx] = numMK * volInv[b];
      numByBox[b] += numMK;
    }

    //Calculate mol fractions
    if (numKinds > 1) {
      for (uint k = 0; k < numKinds; ++k) {
        uint mkIdx = k + numKinds * b;
        if (numByBox[b] > 0) {
          molFractionByKindBox[mkIdx] = (double)numByKindBox[mkIdx] /
                                        (double)numByBox[b];
        } else {
          molFractionByKindBox[mkIdx] = 0.0;
        }
      }
    }

  }

}

#ifdef VARIABLE_PARTICLE_NUMBER


bool MoleculeLookup::ShiftMolBox(const uint mol, const uint currentBox,
                                 const uint intoBox, const uint kind)
{
  uint index = std::find(
                 molLookup + boxAndKindStart[currentBox * numKinds + kind],
                 molLookup + boxAndKindStart[currentBox * numKinds + kind + 1], mol)
               - molLookup;
  assert(index != boxAndKindStart[currentBox * numKinds + kind + 1]);
  assert(molLookup[index] == mol);
  Shift(index, currentBox, intoBox, kind);
  return true;
}

void MoleculeLookup::Shift(const uint index, const uint currentBox,
                           const uint intoBox, const uint kind)
{
  uint oldIndex = index;
  uint newIndex;
  uint section = currentBox * numKinds + kind;
  /* transfer from current box 1 to box 1 */
  /* transfer from current box 1 to box 0 *
  /* transfer from current box 0 to box 0 */  
  boxAndKindSwappableCounts[section]--;
  boxAndKindSwappableCounts[intoBox * numKinds + kind]++;

  if(currentBox >= intoBox) {
    /* Loop invariant - All swappable molecules are 
      first in the section */
    while (section != intoBox * numKinds + kind) {
      newIndex = boxAndKindStart[section]++;
      uint temp = molLookup[oldIndex];
      molLookup[oldIndex] = molLookup[newIndex];
      molLookup[newIndex] = temp;
      /*
      for (int i = 0; i < fixedAtom.size(); i++){
        if (i == NumInBox(0))
          std::cout << " - ";
        //std::cout << fixedAtom[molLookup[i]] << " ";
        std::cout << molLookup[i] << " ";
      }
      std::cout << std::endl;
      */
      SatisfyLoopInvariantDownshift(oldIndex, newIndex, section);
      /*
      for (int i = 0; i < fixedAtom.size(); i++){
        if (i == NumInBox(0))
          std::cout << " - ";
        //std::cout << fixedAtom[molLookup[i]] << " ";
        std::cout << molLookup[i] << " ";
      }
      std::cout << std::endl;
      */
      oldIndex = newIndex;
      --section;
    }
    /* Last iteration - since we add to the back of the section in this if condition */
      /*
      for (int i = 0; i < fixedAtom.size(); i++){
        if (i == NumInBox(0))
          std::cout << " - ";
        //std::cout << fixedAtom[molLookup[i]] << " ";
        std::cout << molLookup[i] << " ";
      }
      std::cout << std::endl;
      */
      SatisfyLoopInvariantDownshift(oldIndex, boxAndKindStart[section]-1, section);
      /*
      for (int i = 0; i < fixedAtom.size(); i++){
        if (i == NumInBox(0))
          std::cout << " - ";
        //std::cout << fixedAtom[molLookup[i]] << " ";
        std::cout << molLookup[i] << " ";
      }
      std::cout << std::endl;
      */
  /* transfer from box 0 to box 1 */  
  /* currentBox < intoBox */
  } else {
    /* Loop invariant - All swappable molecules are 
        first in the section */
    while (section != intoBox * numKinds + kind) {
      newIndex = --boxAndKindStart[++section];
      uint temp = molLookup[oldIndex];
      molLookup[oldIndex] = molLookup[newIndex];
      molLookup[newIndex] = temp;
      /*
      for (int i = 0; i < fixedAtom.size(); i++){
        if (i == NumInBox(0))
          std::cout << " - ";
        //std::cout << fixedAtom[molLookup[i]] << " ";
        std::cout << molLookup[i] << " ";
      }
      std::cout << std::endl;
      */
      SatisfyLoopInvariantUpshift(oldIndex, boxAndKindStart[section-1], section-1);
      /*
      for (int i = 0; i < fixedAtom.size(); i++){
        if (i == NumInBox(0))
          std::cout << " - ";
        //std::cout << fixedAtom[molLookup[i]] << " ";
        std::cout << molLookup[i] << " ";
      }
      std::cout << std::endl;
      */
      oldIndex = newIndex;
    }
    /* No Last Iteration since we add to the front 
       of the section in this else condition */
  }
}

void MoleculeLookup::SatisfyLoopInvariantDownshift(uint oldIndex, uint newIndex, uint section){
  
  /* Edge case - first call : 
      It doesn't enter if condition if the shifted mol was the last swappable or we swapped in place,
      else if there is at least one swappable left it swaps the position of the last swappable 
      (since we decremented boxAndKindSwappableCounts it points to 1 less than the first 
      nonswappable, which is the last swappable).
      with the old Index, which holds the value of the old first value, 
      which by the loop invariant was necessarily swappable.
  */

   /* Edge case - Last call : 
      There is at least 1 swappable since we are adding a molecule, it is in 
      position oldIndex.  We provide boxAndKindStart[section]-1 as the value for newIndex
      (since we incremented boxAndKindSwappableCounts we need to subtract one to get the
      first nonswappable in the section).
  */

   /* This implies a swappable molecule is in the last position and we didn't swap in place*/
  if (boxAndKindSwappableCounts[section] > 0 && oldIndex != newIndex){
    /*
        S* represents the molecule we are shifting down/up
        i.e. boxAndKindSwappableCounts[section] = 2

        BoxAndKindStart
           |
           v
           0      1        2       3   
          swap    swap  nonswap    S* 
        newIndex                oldIndex

      newIndex = boxAndKindStart[section]++;
      uint temp = molLookup[oldIndex];
      molLookup[oldIndex] = molLookup[newIndex];
      molLookup[newIndex] = temp;

              BoxAndKindStart
                  |
                  v
            0     1         2       3   
           S*    swap    nonswap   swap    
        newIndex                 oldIndex

    */

    uint firstNonSwappableIndex = newIndex + boxAndKindSwappableCounts[section];
    uint temp = molLookup[oldIndex];
    molLookup[oldIndex] = molLookup[firstNonSwappableIndex];
    molLookup[firstNonSwappableIndex] = temp;

    /*
                  BoxAndKindStart
                  |
                  v
            0     1     2      3   
           S*    swap  swap  nonswap
        newIndex            oldIndex

    */
  }
}

void MoleculeLookup::SatisfyLoopInvariantUpshift(uint oldIndex, uint newIndex, uint section){
  
  /* Edge case - first call : 
      It doesn't enter if condition if the shifted mol was the last swappable,
      else if there is at least one swappable left it swaps the position of the last swappable 
      (since we decremented boxAndKindSwappableCounts it points to 1 less than the first 
      nonswappable, which is the last swappable).
      with the old Index, which holds the value of the old first value, 
      which by the loop invariant was necessarily swappable.
  */

   /* This implies a swappable molecule is in the last position*/
  if (boxAndKindSwappableCounts[section] > 0){
    /*
        S* represents the molecule we are shifting down/up
        i.e. boxAndKindSwappableCounts[section] = 2

        BoxAndKindStart
           |
           v
           0      1        2       3   
          swap    swap  nonswap    S* 
        newIndex                oldIndex

      newIndex = boxAndKindStart[section]++;
      uint temp = molLookup[oldIndex];
      molLookup[oldIndex] = molLookup[newIndex];
      molLookup[newIndex] = temp;

              BoxAndKindStart
                  |
                  v
            0     1         2       3   
           S*    swap    nonswap   swap    
        newIndex                 oldIndex

    */

    uint firstNonSwappableIndex = newIndex + boxAndKindSwappableCounts[section];
    uint temp = molLookup[oldIndex];
    molLookup[oldIndex] = molLookup[firstNonSwappableIndex];
    molLookup[firstNonSwappableIndex] = temp;

    /*
                  BoxAndKindStart
                  |
                  v
            0     1     2      3   
           S*    swap  swap  nonswap
        newIndex            oldIndex

    */
  }
}


#endif /*ifdef VARIABLE_PARTICLE_NUMBER*/


MoleculeLookup::box_iterator MoleculeLookup::box_iterator::operator++(int)
{
  box_iterator tmp = *this;
  ++(*this);
  return tmp;
}


MoleculeLookup::box_iterator::box_iterator(uint* _pLook, uint* _pSec)
  : pIt(_pLook + * _pSec) {}


MoleculeLookup::box_iterator MoleculeLookup::BoxBegin(const uint box) const
{
  return box_iterator(molLookup, boxAndKindStart + box * numKinds);
}

MoleculeLookup::box_iterator MoleculeLookup::BoxEnd(const uint box) const
{
  return box_iterator(molLookup, boxAndKindStart + (box + 1) * numKinds);
}
