/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Wolf.h"
#include "Ewald.h"
#include "CalculateEnergy.h"
#include "EnergyTypes.h"            //Energy structs
#include "EnsemblePreprocessor.h"   //Flags
#include "BasicTypes.h"             //uint
#include "System.h"                 //For init
#include "StaticVals.h"             //For init
#include "Forcefield.h"
#include "MoleculeKind.h"
#include "Coordinates.h"
#include "BoxDimensions.h"
#include "TrialMol.h"
#include "GeomLib.h"
#include "NumLib.h"
#include <cassert>
#ifdef GOMC_CUDA
#include "CalculateEwaldCUDAKernel.cuh"
#include "CalculateForceCUDAKernel.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#endif
#include "GOMCEventsProfile.h"


Wolf::Wolf(StaticVals & stat, System & sys) :
  Ewald(stat, sys), ffRef(stat.forcefield) {
}

void Wolf::Init() {
    for(uint m = 0; m < mols.count; ++m) {
        const MoleculeKind& molKind = mols.GetKind(m);
        for(uint a = 0; a < molKind.NumAtoms(); ++a) {
          particleKind.push_back(molKind.AtomKind(a));
          particleMol.push_back(m);
          particleCharge.push_back(molKind.AtomCharge(a));
          if(std::abs(molKind.AtomCharge(a)) < 0.000000001) {
              particleHasNoCharge.push_back(true);
          } else {
              particleHasNoCharge.push_back(false);
          }
        }
    }

    // initialize starting index and length index of each molecule
    startMol.resize(currentCoords.Count());
    lengthMol.resize(currentCoords.Count());

    for(int atom = 0; atom < (int) currentCoords.Count(); atom++) {
        startMol[atom] = mols.MolStart(particleMol[atom]);
        lengthMol[atom] = mols.MolLength(particleMol[atom]);
    }

    double molSelfEnergy;
    uint i, j, length;
    molSelfEnergies.resize(mols.kindsCount);
    std::vector<double> & molSelfEnergiesRef = molSelfEnergies;

    for (i = 0; i < mols.GetKindsCount(); i++) {
      MoleculeKind const& thisKind = mols.kinds[i];
      length = thisKind.NumAtoms();
      molSelfEnergy = 0.0;
      for (j = 0; j < length; j++) {
          molSelfEnergy += (thisKind.AtomCharge(j) * thisKind.AtomCharge(j));
      }
      // Store this for quick access in ChangeSelf
      molSelfEnergiesRef[i] = molSelfEnergy;
    }

    oneThree = ff.OneThree;
    oneFour = ff.OneFour;
    scaling_14 = ff.scaling_14;
    SetCoulKind(ff.coulKind);
    SetWolfKind(ff.wolfKind);
}

void Wolf::SetWolfKind(uint wolfKindArg){
    wolfKind = wolfKindArg;
    switch(wolfKindArg) {
      //WOLF_HYBRID_KIND
      case 0:
        isVlugtWolf = false;
        isGrossWolf = false;
        isHybridWolf = true;
        isVlugtWithIntraCutoffWolf = false;
        break;
      // WOLF_VLUGT_KIND
      case 1:
        isVlugtWolf = true;
        isGrossWolf = false;
        isHybridWolf = false;
        isVlugtWithIntraCutoffWolf = false;
        break;
      // WOLF_GROSS_KIND
      case 2:
        isVlugtWolf = false;
        isGrossWolf = true;
        isHybridWolf = false;
        isVlugtWithIntraCutoffWolf = false;
        break;
      case 3:
        isVlugtWolf = false;
        isGrossWolf = false;
        isHybridWolf = false;
        isVlugtWithIntraCutoffWolf = true;
        break;
      default:
        std::cout << "Error ff.WolfKind has invalid value!  Check WolfKind in Config File!" << std::endl;
        exit(1);
    }
}

void Wolf::SetCoulKind(uint coulKindArg){
  if (coulKindArg > COUL_TOTAL_KINDS){
    std::cout << "Error ff.coulKind has invalid value!  Check coulKind in Config File!" << std::endl;
    exit(1);
  } else{
    coulKind = coulKindArg;
  }
}

uint Wolf::GetWolfKind(void){
  return wolfKind;
}

uint Wolf::GetCoulKind(void){
  return coulKind;
}



void Wolf::AllocMem()
{
  return;
}

void Wolf::RecipInit(uint box, BoxDimensions const& boxAxes)
{
  return;
}

//compute reciprocal term for a box with a new volume
void Wolf::BoxReciprocalSetup(uint box, XYZArray const& molCoords)
{
  return;
}

//compute reciprocal term for a box when not testing a volume change
void Wolf::BoxReciprocalSums(uint box, XYZArray const& molCoords)
{
  return;
}

//calculate reciprocal term for a box
double Wolf::BoxReciprocal(uint box, bool isNewVolume) const
{
  return 0.0;
}

//calculate reciprocal force term for a box with molCoords
void Wolf::BoxForceReciprocal(XYZArray const& molCoords,
                                 XYZArray& atomForceRec, XYZArray& molForceRec,
                                 uint box)
{
  return;
}

//calculate reciprocal force term for a box
Virial Wolf::VirialReciprocal(Virial& virial, uint box) const
{
  return virial;
}


//calculate reciprocal term for displacement and rotation move
double Wolf::MolReciprocal(XYZArray const& molCoords,
                              const uint molIndex, const uint box)
{
  return 0.0;
}


//calculate reciprocal term for lambdaNew and Old with same coordinates
double Wolf::ChangeLambdaRecip(XYZArray const& molCoords, const double lambdaOld,
                                  const double lambdaNew, const uint molIndex,
                                  const uint box)
{
  return 0.0;
}

//calculate self term for a box
double Wolf::BoxSelf(uint box) const
{
  return BoxSelf(box, ffRef.wolfFactor1[box], ffRef.wolfAlpha[box]);
}

//calculate self term for a box
double Wolf::BoxSelf(uint box,
                      double wolfFactor1,
                      double wolfAlpha) const
{
    if (box >= BOXES_WITH_U_NB){
        return 0.0;
    } else {
        GOMC_EVENT_START(1, GomcProfileEvent::SELF_BOX);
        double self = 0.0;
        double molSelfEnergy;
        uint i, j, length, molNum;
        double lambdaCoef = 1.0;
        const std::vector<double> & molSelfEnergiesRef = molSelfEnergies;

        for (i = 0; i < mols.GetKindsCount(); i++) {
          MoleculeKind const& thisKind = mols.kinds[i];
          molNum = molLookup.NumKindInBox(i, box);
          molSelfEnergy = 0.0;
          if(lambdaRef.KindIsFractional(i, box)) {
              //If a molecule is fractional, we subtract the fractional molecule and
              // add it later
              --molNum;
              //returns lambda and not sqrt(lambda)
              lambdaCoef = lambdaRef.GetLambdaCoulomb(i, box);
          }
          // Store this for quick access in ChangeSelf
          molSelfEnergy = molSelfEnergiesRef[i];
          self += (molSelfEnergy * molNum);
          if(lambdaRef.KindIsFractional(i, box)) {
              //Add the fractional molecule part
              self += (molSelfEnergy * lambdaCoef);
          }
        }
        if (isVlugtWolf || isVlugtWithIntraCutoffWolf){
          self *= -0.5 * ((wolfAlpha * M_2_SQRTPI) + wolfFactor1);
        } else {
          // we eliminate the alpha/root(pi) using Wolf,mod from Gross et al
          self *= -0.5 * wolfFactor1;
        }

        GOMC_EVENT_STOP(1, GomcProfileEvent::SELF_BOX);
        return self * num::qqFact;
    }
}

//calculate self term for a box
double Wolf::MolCorrection(uint molIndex, uint box) const
{
  return MolCorrection(molIndex, box, ffRef.rCutCoulomb[box], ffRef.rCutCoulombSq[box],
                                      ffRef.wolfFactor1[box], 
                                      ffRef.wolfFactor2[box], 
                                      ffRef.wolfAlpha[box]);
}

//calculate correction term for a molecule
double Wolf::MolCorrection(uint molIndex, uint box,
                          double rCutCoulomb,
                          double rCutCoulombSq,
                          double wolfFactor1,
                          double wolfFactor2,
                          double wolfAlpha) const
{
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::CORR_MOL);
  double dist, distSq;
  double correction = 0.0, dampenedCorr = 0.0, undampenedCorr = 0.0;
  XYZ virComponents;

  MoleculeKind& thisKind = mols.kinds[mols.kIndex[molIndex]];
  uint atomSize = thisKind.NumAtoms();
  uint start = mols.MolStart(molIndex);
  double lambdaCoef = GetLambdaCoef(molIndex, box);


  // This term only needs to be calculated once, if the molecules are rigid
  // Or of the maximum radius of gyration < RCutCoulomb
  // Otherwise, parts of the molecule may extend out of range of each other
  // For now assume the latter.
  for (uint i = 0; i < atomSize; i++) {
    if(particleHasNoCharge[start + i]) {
      continue;
    }
    for (uint j = i + 1; j < atomSize; j++) {
      currentAxes.InRcut(distSq, virComponents, currentCoords,
                        start + i, start + j, box);
      dampenedCorr = 0.0;
      if(distSq < rCutCoulombSq){
        if (!isHybridWolf)
          dampenedCorr -= wolfFactor1;
        if (isVlugtWolf || isVlugtWithIntraCutoffWolf){
          dist = sqrt(distSq);
          dampenedCorr += -1.0*erf(wolfAlpha * dist)/dist;   
          if(ff.coulKind && isVlugtWithIntraCutoffWolf){
            double distDiff = dist-rCutCoulomb;
            dampenedCorr += wolfFactor2*distDiff;
          } 
        }
        correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * dampenedCorr;
      }
    }
  }
  if (isGrossWolf || isHybridWolf){
    for (uint i = 0; i < atomSize; i++) {
      if(oneThree){
        //loop over all 1-3 partners of the particle
        const uint* partner = thisKind.sortedNB_1_3.Begin(i);
        const uint* end = thisKind.sortedNB_1_3.End(i);
        while (partner != end) {
          // Need to check for cutoff for all kinds
          currentAxes.InRcut(distSq, virComponents, currentCoords,
                            start + i, start + (*partner), box); 
          if(distSq < rCutCoulombSq){
            dist = sqrt(distSq);
            if (isGrossWolf){
              // scale the erfc term
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr *= scaling_14;
            } else if (isHybridWolf) {
              // scale the entire erfc term and wolf factor 1
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr -= wolfFactor1;
              dampenedCorr *= scaling_14;
            }
            correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
          }
          ++partner;
        }
      }
      if(oneFour){
        //loop over all 1-4 partners of the particle
        const uint* partner = thisKind.sortedNB_1_4.Begin(i);
        const uint* end = thisKind.sortedNB_1_4.End(i);
        while (partner != end) {
          // Need to check for cutoff for all kinds
          currentAxes.InRcut(distSq, virComponents, currentCoords,
                            start + i, start + (*partner), box); 
          if(distSq < rCutCoulombSq){
            dist = sqrt(distSq);
            if (isGrossWolf){
              // scale the erfc term
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr *= scaling_14;
            } else if (isHybridWolf) {
              // scale the entire erfc term and wolf factor 1
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr -= wolfFactor1;
              dampenedCorr *= scaling_14;
            }
            correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
          }
          ++partner;
        }
      }
      //loop over all 1-N partners of the particle
      const uint* partner = thisKind.sortedNB.Begin(i);
      const uint* end = thisKind.sortedNB.End(i);
      while (partner != end) {
        // Need to check for cutoff for all kinds
        currentAxes.InRcut(distSq, virComponents, currentCoords,
                            start + i, start + (*partner), box); 
        if(distSq < rCutCoulombSq){
          dist = sqrt(distSq);
          if (isGrossWolf){
            // scale the erfc term
            dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist; 
            // Don't scale 1-N 
            //dampenedCorr *= scaling_14;
          } else if (isHybridWolf) {
            // scale the entire erfc term and wolf factor 1
            dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
            dampenedCorr -= wolfFactor1;
            // Don't scale 1-N
            //dampenedCorr *= scaling_14;
          }
          correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
        }
        ++partner;
      }      
    } 
  }
  

  GOMC_EVENT_STOP(1, GomcProfileEvent::CORR_MOL);
  return num::qqFact * correction * lambdaCoef * lambdaCoef;
}

//calculate reciprocal term in destination box for swap move
double Wolf::SwapDestRecip(const cbmc::TrialMol &newMol,
                              const uint box,
                              const int molIndex)
{
  return 0.0;
}


//calculate reciprocal term in source box for swap move
double Wolf::SwapSourceRecip(const cbmc::TrialMol &oldMol,
                                const uint box, const int molIndex)
{
  return 0.0;
}


//calculate reciprocal term for inserting some molecules (kindA) in destination
// box and removing a molecule (kindB) from destination box
double Wolf::MolExchangeReciprocal(const std::vector<cbmc::TrialMol> &newMol,
                                      const std::vector<cbmc::TrialMol> &oldMol,
                                      const std::vector<uint> &molIndexNew,
                                      const std::vector<uint> &molIndexold,
                                      bool first_call)
{
  return 0.0;
}

double Wolf::SwapSelf(const cbmc::TrialMol& trialMol) const{
  uint box = trialMol.GetBox();
  return SwapSelf(trialMol,
                      ffRef.wolfFactor1[box],
                      ffRef.wolfAlpha[box]);
}

//calculate self term after swap move
double Wolf::SwapSelf(const cbmc::TrialMol& trialMol,
                      double wolfFactor1,
                      double wolfAlpha) const
{
  uint box = trialMol.GetBox();
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::SELF_SWAP);
  MoleculeKind const& thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();
  double en_self = 0.0;

  en_self = molSelfEnergies[thisKind.kindIndex];

  GOMC_EVENT_STOP(1, GomcProfileEvent::SELF_SWAP);
  if (isVlugtWolf || isVlugtWithIntraCutoffWolf){
    en_self *= (((wolfAlpha * M_2_SQRTPI) + wolfFactor1) * -0.5);
  } else {
    // we eliminate the alpha/root(pi) using Wolf,mod from Gross et al
    en_self *= (wolfFactor1 * -0.5);
  }
  return en_self  * num::qqFact;
}


//calculate correction term for a molecule with system lambda
//It's called when the molecule configuration changes, regrowth, crankshaft, IntraSwap, IntraMEMC ...
double Wolf::SwapCorrection(const cbmc::TrialMol& trialMol) const {
  uint box = trialMol.GetBox();
  return SwapCorrection(trialMol,
                        box, ffRef.rCutCoulomb[box], 
                        ffRef.rCutCoulombSq[box],
                        ffRef.wolfFactor1[box], 
                        ffRef.wolfFactor2[box], 
                        ffRef.wolfAlpha[box]);
}

//calculate correction term for a molecule with system lambda
//It's called when the molecule configuration changes, regrowth, crankshaft, IntraSwap, IntraMEMC ...
double Wolf::SwapCorrection(const cbmc::TrialMol& trialMol,
                             const uint molIndex) const 
{
  uint box = trialMol.GetBox();
  return SwapCorrection(trialMol,
                        molIndex, 
                        box, ffRef.rCutCoulomb[box], 
                        ffRef.rCutCoulombSq[box],
                        ffRef.wolfFactor1[box], 
                        ffRef.wolfFactor2[box], 
                        ffRef.wolfAlpha[box]);
}

//calculate correction term after swap move
double Wolf::SwapCorrection(const cbmc::TrialMol& trialMol,
                          const uint molIndex,
                          uint box,
                          double rCutCoulomb,
                          double rCutCoulombSq,
                          double wolfFactor1,
                          double wolfFactor2,
                          double wolfAlpha) const
{
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::CORR_SWAP);
  double dist, distSq;
  double correction = 0.0, dampenedCorr = 0.0, undampenedCorr = 0.0;
  XYZ virComponents;

  MoleculeKind& thisKind = mols.kinds[mols.kIndex[molIndex]];
  uint atomSize = thisKind.NumAtoms();
  uint start = mols.MolStart(molIndex);
  double lambdaCoef = GetLambdaCoef(molIndex, box);


    // This term only needs to be calculated once, if the molecules are rigid
    // Or of the maximum radius of gyration < RCutCoulomb
    // Otherwise, parts of the molecule may extend out of range of each other
    // For now assume the latter.
  for (uint i = 0; i < atomSize; i++) {
    if(particleHasNoCharge[start + i]) {
      continue;
    }
    for (uint j = i + 1; j < atomSize; j++) {
      currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                         i, j, box);
      dampenedCorr = 0.0;
      if(distSq < rCutCoulombSq){
        if (!isHybridWolf)
          dampenedCorr -= wolfFactor1;
        if (isVlugtWolf || isVlugtWithIntraCutoffWolf){
          dist = sqrt(distSq);
          dampenedCorr += -1.0*erf(wolfAlpha * dist)/dist;   
          if(ff.coulKind && isVlugtWithIntraCutoffWolf){
            double distDiff = dist-rCutCoulomb;
            dampenedCorr += wolfFactor2*distDiff;
          } 
        }
        correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * dampenedCorr;
      }
    }
  }
  if (isGrossWolf || isHybridWolf){
    for (uint i = 0; i < atomSize; i++) {
      if(oneThree){
        //loop over all 1-3 partners of the particle
        const uint* partner = thisKind.sortedNB_1_3.Begin(i);
        const uint* end = thisKind.sortedNB_1_3.End(i);
        while (partner != end) {
          // Need to check for cutoff for all kinds
          currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                         i, *partner, box);
          if(distSq < rCutCoulombSq){
            dist = sqrt(distSq);
            if (isGrossWolf){
              // scale the erfc term
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr *= scaling_14;
            } else if (isHybridWolf) {
              // scale the entire erfc term and wolf factor 1
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr -= wolfFactor1;
              dampenedCorr *= scaling_14;
            }
            correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
          }
          ++partner;
        }
      }
      if(oneFour){
        //loop over all 1-4 partners of the particle
        const uint* partner = thisKind.sortedNB_1_4.Begin(i);
        const uint* end = thisKind.sortedNB_1_4.End(i);
        while (partner != end) {
          // Need to check for cutoff for all kinds
          currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                         i, *partner, box);
          if(distSq < rCutCoulombSq){
            dist = sqrt(distSq);
            if (isGrossWolf){
              // scale the erfc term
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr *= scaling_14;
            } else if (isHybridWolf) {
              // scale the entire erfc term and wolf factor 1
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr -= wolfFactor1;
              dampenedCorr *= scaling_14;
            }
            correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
          }
          ++partner;
        }
      }
      //loop over all 1-N partners of the particle
      const uint* partner = thisKind.sortedNB.Begin(i);
      const uint* end = thisKind.sortedNB.End(i);
      while (partner != end) {
        // Need to check for cutoff for all kinds
        currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                        i, *partner, box);
        if(distSq < rCutCoulombSq){
          dist = sqrt(distSq);
          if (isGrossWolf){
            // scale the erfc term
            dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist; 
            // Don't scale 1-N 
            //dampenedCorr *= scaling_14;
          } else if (isHybridWolf) {
            // scale the entire erfc term and wolf factor 1
            dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
            dampenedCorr -= wolfFactor1;
            // Don't scale 1-N
            //dampenedCorr *= scaling_14;
          }
          correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
        }
        ++partner;
      }      
    } 
  }
  

  GOMC_EVENT_STOP(1, GomcProfileEvent::CORR_SWAP);
  return num::qqFact * correction * lambdaCoef * lambdaCoef;
}


//calculate correction term after swap move
double Wolf::SwapCorrection(const cbmc::TrialMol& trialMol,
                          uint box,
                          double rCutCoulomb,
                          double rCutCoulombSq,
                          double wolfFactor1,
                          double wolfFactor2,
                          double wolfAlpha) const
{
  if (box >= BOXES_WITH_U_NB)
    return 0.0;

  GOMC_EVENT_START(1, GomcProfileEvent::CORR_SWAP);
  double dist, distSq;
  double correction = 0.0, dampenedCorr = 0.0, undampenedCorr = 0.0;
  XYZ virComponents;

  const MoleculeKind& thisKind = trialMol.GetKind();
  uint atomSize = thisKind.NumAtoms();


    // This term only needs to be calculated once, if the molecules are rigid
    // Or of the maximum radius of gyration < RCutCoulomb
    // Otherwise, parts of the molecule may extend out of range of each other
    // For now assume the latter.
  for (uint i = 0; i < atomSize; i++) {
    for (uint j = i + 1; j < atomSize; j++) {
      currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                         i, j, box);
      dampenedCorr = 0.0;
      if(distSq < rCutCoulombSq){
        if (!isHybridWolf)
          dampenedCorr -= wolfFactor1;
        if (isVlugtWolf || isVlugtWithIntraCutoffWolf){
          dist = sqrt(distSq);
          dampenedCorr += -1.0*erf(wolfAlpha * dist)/dist;   
          if(ff.coulKind && isVlugtWithIntraCutoffWolf){
            double distDiff = dist-rCutCoulomb;
            dampenedCorr += wolfFactor2*distDiff;
          } 
        }
        correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * dampenedCorr;
      }
    }
  }
  if (isGrossWolf || isHybridWolf){
    for (uint i = 0; i < atomSize; i++) {
      if(oneThree){
        //loop over all 1-3 partners of the particle
        const uint* partner = thisKind.sortedNB_1_3.Begin(i);
        const uint* end = thisKind.sortedNB_1_3.End(i);
        while (partner != end) {
          // Need to check for cutoff for all kinds
          currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                         i, *partner, box); 
          if(distSq < rCutCoulombSq){
            dist = sqrt(distSq);
            if (isGrossWolf){
              // scale the erfc term
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr *= scaling_14;
            } else if (isHybridWolf) {
              // scale the entire erfc term and wolf factor 1
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr -= wolfFactor1;
              dampenedCorr *= scaling_14;
            }
            correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
          }
          ++partner;
        }
      }
      if(oneFour){
        //loop over all 1-4 partners of the particle
        const uint* partner = thisKind.sortedNB_1_4.Begin(i);
        const uint* end = thisKind.sortedNB_1_4.End(i);
        while (partner != end) {
          // Need to check for cutoff for all kinds
          currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                         i, *partner, box); 
          if(distSq < rCutCoulombSq){
            dist = sqrt(distSq);
            if (isGrossWolf){
              // scale the erfc term
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr *= scaling_14;
            } else if (isHybridWolf) {
              // scale the entire erfc term and wolf factor 1
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr -= wolfFactor1;
              dampenedCorr *= scaling_14;
            }
            correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
          }
          ++partner;
        }
      }
      //loop over all 1-N partners of the particle
      const uint* partner = thisKind.sortedNB.Begin(i);
      const uint* end = thisKind.sortedNB.End(i);
      while (partner != end) {
        // Need to check for cutoff for all kinds
        currentAxes.InRcut(distSq, virComponents, trialMol.GetCoords(),
                         i, *partner, box); 
        if(distSq < rCutCoulombSq){
          dist = sqrt(distSq);
          if (isGrossWolf){
            // scale the erfc term
            dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist; 
            // Don't scale 1-N 
            //dampenedCorr *= scaling_14;
          } else if (isHybridWolf) {
            // scale the entire erfc term and wolf factor 1
            dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
            dampenedCorr -= wolfFactor1;
            // Don't scale 1-N
            //dampenedCorr *= scaling_14;
          }
          correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
        }
        ++partner;
      }      
    } 
  }
  

  GOMC_EVENT_STOP(1, GomcProfileEvent::CORR_SWAP);
  return num::qqFact * correction;
}


//It's called in free energy calculation to calculate the change in
// self energy in all lambda states
void Wolf::ChangeSelf(Energy *energyDiff, Energy &dUdL_Coul,
                          const std::vector<double> &lambda_Coul,
                          const uint iState, const uint molIndex,
                          const uint box) const{
  ChangeSelf(energyDiff, dUdL_Coul, lambda_Coul, iState, molIndex,
              box, ffRef.wolfFactor1[box], ffRef.wolfAlpha[box]);
}

//It's called in free energy calculation to calculate the change in
// self energy in all lambda states
void Wolf::ChangeSelf(Energy *energyDiff, Energy &dUdL_Coul,
                         const std::vector<double> &lambda_Coul,
                         const uint iState, const uint molIndex,
                         const uint box,
                          double wolfFactor1,
                          double wolfAlpha) const
{
    //Vlugt
    //en_self *= ((wolfAlpha * M_2_SQRTPI) +  wolfFactor1 );
    // We eliminate the alpha/root(pi) using Wolf,mod

    uint atomSize = mols.GetKind(molIndex).NumAtoms();
    uint start = mols.MolStart(molIndex);
    uint lambdaSize = lambda_Coul.size();
    double coefDiff, en_self = 0.0;
    //Calculate the self energy with lambda = 1
    for (uint i = 0; i < atomSize; i++) {
      en_self += (particleCharge[i + start] * particleCharge[i + start]);
    }
    if (isVlugtWolf || isVlugtWithIntraCutoffWolf){
      en_self *= -0.5 * ((wolfAlpha * M_2_SQRTPI) + wolfFactor1) * num::qqFact;
    } else {
      // we eliminate the alpha/root(pi) using Wolf,mod from Gross et al
      en_self *= -0.5 * wolfFactor1 * num::qqFact;
    }

    //Calculate the energy difference for each lambda state
    for (uint s = 0; s < lambdaSize; s++) {
      coefDiff = lambda_Coul[s] - lambda_Coul[iState];
      energyDiff[s].self += coefDiff * en_self;
    }
    //Calculate du/dl of self for current state, for linear scaling
    dUdL_Coul.self += en_self;
}

//It's called in free energy calculation to calculate the change in
// correction energy in all lambda states
void Wolf::ChangeCorrection(Energy *energyDiff, Energy &dUdL_Coul,
                             const std::vector<double> &lambda_Coul,
                             const uint iState, const uint molIndex,
                             const uint box) const{
  ChangeCorrection(energyDiff, dUdL_Coul, lambda_Coul, iState,
                  molIndex, box, ffRef.rCutCoulomb[box], ffRef.rCutCoulombSq[box],
                                      ffRef.wolfFactor1[box], 
                                      ffRef.wolfFactor2[box], 
                                      ffRef.wolfAlpha[box]);
}

//It's called in free energy calculation to calculate the change in
// correction energy in all lambda states
void Wolf::ChangeCorrection(Energy *energyDiff, Energy &dUdL_Coul,
                               const std::vector<double> &lambda_Coul,
                               const uint iState, const uint molIndex,
                               const uint box,
                                double rCutCoulomb,
                                double rCutCoulombSq,
                                double wolfFactor1,
                                double wolfFactor2,
                                double wolfAlpha) const
{
  uint atomSize = mols.GetKind(molIndex).NumAtoms();
  uint start = mols.MolStart(molIndex);
  uint lambdaSize = lambda_Coul.size();
  double coefDiff, distSq, dist, correction = 0.0, dampenedCorr;
  XYZ virComponents;
  MoleculeKind& thisKind = mols.kinds[mols.kIndex[molIndex]];


  //Calculate the correction energy with lambda = 1
  for (uint i = 0; i < atomSize; i++) {
    if(particleHasNoCharge[start + i]) {
      continue;
    }
    for (uint j = i + 1; j < atomSize; j++) {
      currentAxes.InRcut(distSq, virComponents, currentCoords,
                        start + i, start + j, box);
      dampenedCorr = 0.0;
      if(distSq < rCutCoulombSq){
        if (!isHybridWolf)
          dampenedCorr -= wolfFactor1;
        if (isVlugtWolf || isVlugtWithIntraCutoffWolf){
          dist = sqrt(distSq);
          dampenedCorr += -1.0*erf(wolfAlpha * dist)/dist;   
          if(ff.coulKind && isVlugtWithIntraCutoffWolf){
            double distDiff = dist-rCutCoulomb;
            dampenedCorr += wolfFactor2*distDiff;
          } 
        }
        correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(j) * dampenedCorr;
      }
    }
  }
  if (isGrossWolf || isHybridWolf){
    for (uint i = 0; i < atomSize; i++) {
      if(oneThree){
        //loop over all 1-3 partners of the particle
        const uint* partner = thisKind.sortedNB_1_3.Begin(i);
        const uint* end = thisKind.sortedNB_1_3.End(i);
        while (partner != end) {
          // Need to check for cutoff for all kinds
          currentAxes.InRcut(distSq, virComponents, currentCoords,
                            start + i, start + (*partner), box); 
          if(distSq < rCutCoulombSq){
            dist = sqrt(distSq);
            if (isGrossWolf){
              // scale the erfc term
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr *= scaling_14;
            } else if (isHybridWolf) {
              // scale the entire erfc term and wolf factor 1
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr -= wolfFactor1;
              dampenedCorr *= scaling_14;
            }
            correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
          }
          ++partner;
        }
      }
      if(oneFour){
        //loop over all 1-4 partners of the particle
        const uint* partner = thisKind.sortedNB_1_4.Begin(i);
        const uint* end = thisKind.sortedNB_1_4.End(i);
        while (partner != end) {
          // Need to check for cutoff for all kinds
          currentAxes.InRcut(distSq, virComponents, currentCoords,
                            start + i, start + (*partner), box); 
          if(distSq < rCutCoulombSq){
            dist = sqrt(distSq);
            if (isGrossWolf){
              // scale the erfc term
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr *= scaling_14;
            } else if (isHybridWolf) {
              // scale the entire erfc term and wolf factor 1
              dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
              dampenedCorr -= wolfFactor1;
              dampenedCorr *= scaling_14;
            }
            correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
          }
          ++partner;
        }
      }
      //loop over all 1-N partners of the particle
      const uint* partner = thisKind.sortedNB.Begin(i);
      const uint* end = thisKind.sortedNB.End(i);
      while (partner != end) {
        // Need to check for cutoff for all kinds
        currentAxes.InRcut(distSq, virComponents, currentCoords,
                          start + i, start + (*partner), box); 
        if(distSq < rCutCoulombSq){
          dist = sqrt(distSq);
          if (isGrossWolf){
            // scale the erfc term
            dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist; 
            // Don't scale 1-N 
            //dampenedCorr *= scaling_14;
          } else if (isHybridWolf) {
            // scale the entire erfc term and wolf factor 1
            dampenedCorr = -1.0*erf(wolfAlpha * dist)/dist;  
            dampenedCorr -= wolfFactor1;
            // Don't scale 1-N
            //dampenedCorr *= scaling_14;
          }
          correction += thisKind.AtomCharge(i) * thisKind.AtomCharge(*partner) * dampenedCorr;
        }
        ++partner;
      }      
    } 
  }
  
  //Calculate the energy difference for each lambda state
  for (uint s = 0; s < lambdaSize; s++) {
    coefDiff = lambda_Coul[s] - lambda_Coul[iState];
    energyDiff[s].correction += coefDiff * correction;
  }
  //Calculate du/dl of correction for current state, for linear scaling
  dUdL_Coul.correction += correction;
}

//It's called in free energy calculation to calculate the change in
// reciprocal energy in all lambda states
void Wolf::ChangeRecip(Energy *energyDiff, Energy &dUdL_Coul,
                          const std::vector<double> &lambda_Coul,
                          const uint iState, const uint molIndex,
                          const uint box) const
{
  return;
}


//back up reciprocal value to Ref (will be called during initialization)
void Wolf::SetRecipRef(uint box)
{
  return;
}

//update reciprocal values
void Wolf::UpdateRecip(uint box)
{
  return;
}

//copy reciprocal values from ref to new
void Wolf::CopyRecip(uint box)
{
  return;
}

//update the kx, ky, kz, hsqr and prefact
void Wolf::UpdateRecipVec(uint box)
{
  return;
}

//restore cosMol and sinMol
void Wolf::RestoreMol(int molIndex)
{
  return;
}

//restore the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void Wolf::exgMolCache()
{
  return;
}

//backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void Wolf::backupMolCache()
{
  return;
}

void Wolf::UpdateVectorsAndRecipTerms(bool output)
{
  return;
}
