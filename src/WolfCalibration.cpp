#include "WolfCalibration.h"
WolfCalibration::WolfCalibration(config_setup::WolfCalibration const& wolfCal,
                                  bool enable)
{
  if(enable){
    CalculateWolfCalibrationMemoryUsage(wolfCal);
    AllocMem();
    InitWolfCalibration(wolfCal);
  }
}

WolfCalibration::~WolfCalibration()
{
  DeallocMem();
}

double WolfCalibration::GetAlpha(int box, int indexForAlpha){return wolfAlpha[startOfNumAlphas[box]+indexForAlpha];}
double WolfCalibration::GetRCut(int box, int indexForRCut){return rCutCoulomb[startOfNumRCuts[box]+indexForRCut];}
double WolfCalibration::GetRCutSq(int box, int indexForRCut){return rCutCoulombSq[startOfNumRCuts[box]+indexForRCut];}
double WolfCalibration::GetWolfFactor1(int box, int indexForRCut, int indexForAlpha){return wolfFactor1[startOfWolfFactor[box] + numberOfRCuts[box]*indexForRCut + indexForAlpha];}
double WolfCalibration::GetWolfFactor2(int box, int indexForRCut, int indexForAlpha){return wolfFactor2[startOfWolfFactor[box] + numberOfRCuts[box]*indexForRCut + indexForAlpha];}
double WolfCalibration::GetWolfFactor3(int box, int indexForRCut, int indexForAlpha){return wolfFactor3[startOfWolfFactor[box] + numberOfRCuts[box]*indexForRCut + indexForAlpha];}
double WolfCalibration::GetNumberOfRCuts(int box){return numberOfRCuts[box];}
double WolfCalibration::GetNumberOfAlphas(int box){return numberOfAlphas[box];}
int WolfCalibration::GetTotalNumWolfFactors(){return totalNumWolfFactors;}
int WolfCalibration::GetStartOfWolfFactors(int box){return startOfWolfFactor[box];}

void WolfCalibration::InitWolfCalibration(config_setup::WolfCalibration const& wolfCal){
  for(uint b = 0 ; b < BOX_TOTAL; b++) {
    for(uint r = 0; r < numberOfRCuts[b]; r++) {
      rCutCoulomb[startOfNumRCuts[b]+r] = wolfCal.wolfCutoffCoulombStart[b] + r*wolfCal.wolfCutoffCoulombDelta[b];
      rCutCoulombSq[startOfNumRCuts[b]+r] = rCutCoulomb[startOfNumRCuts[b]+r] * rCutCoulomb[startOfNumRCuts[b]+r];
    }
    for(uint a = 0 ; a < numberOfAlphas[b]; a++) {
      wolfAlpha[startOfNumAlphas[b]+a] = wolfCal.wolfAlphaStart[b] + a*wolfCal.wolfAlphaDelta[b];
    }
    for(uint r = 0; r < numberOfRCuts[b]; r++) {
      for(uint a = 0 ; a < numberOfAlphas[b]; a++) {
        wolfFactor1[startOfWolfFactor[b] + numberOfRCuts[b]*r + a] = erfc(wolfAlpha[startOfNumAlphas[b]+a]*rCutCoulomb[startOfNumRCuts[b]+r])/rCutCoulomb[startOfNumRCuts[b]+r];
        wolfFactor2[startOfWolfFactor[b] + numberOfRCuts[b]*r + a] = wolfFactor1[startOfWolfFactor[b] + numberOfRCuts[b]*r + a]/rCutCoulomb[startOfNumRCuts[b]+r];
        wolfFactor2[startOfWolfFactor[b] + numberOfRCuts[b]*r + a] += wolfAlpha[startOfNumAlphas[b]+a] *  M_2_SQRTPI * 
                          exp(-1.0*wolfAlpha[startOfNumAlphas[b]+a]*wolfAlpha[startOfNumAlphas[b]+a]*rCutCoulombSq[startOfNumRCuts[b]+r])
                          /rCutCoulomb[startOfNumRCuts[b]+r];
        wolfFactor3[startOfWolfFactor[b] + numberOfRCuts[b]*r + a] = wolfAlpha[startOfNumAlphas[b]+a] *  M_2_SQRTPI;
      }
    }
  }
}

void WolfCalibration::CalculateWolfCalibrationMemoryUsage(config_setup::WolfCalibration const& wolfCal){
  // Calculate the number of calibration points from the ranges provided
  // If delta is not a common multiple of the (End - Start) explicitly add the End.
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    numberOfRCuts[b] = (int)((wolfCal.wolfCutoffCoulombEnd[b] - wolfCal.wolfCutoffCoulombStart[b]) / wolfCal.wolfCutoffCoulombDelta[b]) + 1;
    numberOfAlphas[b] = (int)((wolfCal.wolfAlphaEnd[b] - wolfCal.wolfAlphaStart[b]) / wolfCal.wolfAlphaDelta[b]) + 1;
  }
}

void WolfCalibration::AllocMem(){
  totalNumAlphas = 0;
  totalNumRCuts = 0;
  totalNumWolfFactors = 0;
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    startOfNumAlphas[b] = totalNumAlphas;
    startOfNumRCuts[b] = totalNumRCuts;
    startOfWolfFactor[b] = totalNumWolfFactors;
    totalNumAlphas += numberOfAlphas[b];
    totalNumRCuts += numberOfRCuts[b];
    totalNumWolfFactors +=  numberOfAlphas[b] * numberOfRCuts[b];
  }

  wolfAlpha = new double[totalNumAlphas];
  rCutCoulomb = new double[totalNumRCuts];
  rCutCoulombSq = new double[totalNumRCuts];
  wolfFactor1 = new double[totalNumWolfFactors];
  wolfFactor2 = new double[totalNumWolfFactors];
  wolfFactor3 = new double[totalNumWolfFactors];
}

void WolfCalibration::DeallocMem(){
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    delete[] wolfAlpha;
    delete[] rCutCoulomb;
    delete[] rCutCoulombSq;
    delete[] wolfFactor1;
    delete[] wolfFactor2;
    delete[] wolfFactor3;
  }
}

int WolfCalibration::GetIndex(int box, int wolfKind, int coulKind, int r, int a){
                  // Points to start of array if box 0
                  // Skips box 0 entries if box 1, 
      int index = WOLF_TOTAL_KINDS*COUL_TOTAL_KINDS*startOfWolfFactor[box]
                  // Skip wolf kind entries
                  + wolfKind*COUL_TOTAL_KINDS*numberOfRCuts[box]*numberOfAlphas[box]
                  // skip coul kind entries
                  + coulKind*numberOfRCuts[box]*numberOfAlphas[box]
                  // Within a file (WOLFKIND_COULKIND_BOX)
                  + startOfWolfFactor[box] + r*numberOfAlphas[box] + a;
      return index;
}