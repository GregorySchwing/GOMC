#include "WolfCalibration.h"
WolfCalibration::WolfCalibration(config_setup::WolfCalibration const& wolfCal){
  CalculateWolfCalibrationMemoryUsage(wolfCal);
  AllocMem();
  InitWolfCalibration(wolfCal);
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

void WolfCalibration::InitWolfCalibration(config_setup::WolfCalibration const& wolfCal){
  for(uint b = 0 ; b < BOX_TOTAL; b++) {
    for(uint r = 0; r < numberOfRCuts[b]-1; r++) {
      rCutCoulomb[startOfNumRCuts[b]+r] = wolfCal.wolfCutoffCoulombStart[b] + r*wolfCal.wolfCutoffCoulombDelta[b];
      rCutCoulombSq[startOfNumRCuts[b]+r] = rCutCoulomb[startOfNumRCuts[b]+r] * rCutCoulomb[startOfNumRCuts[b]+r];
    }
    for(uint a = 0 ; a < numberOfAlphas[b]-1; a++) {
      wolfAlpha[startOfNumAlphas[b]+a] = wolfCal.wolfAlphaStart[b] + a*wolfCal.wolfAlphaDelta[b];
    }
    if (explicitlyAddEndRCut[b]){
      rCutCoulomb[startOfNumRCuts[b]+numberOfRCuts[b]-1] = wolfCal.wolfCutoffCoulombEnd[b];
      rCutCoulombSq[startOfNumRCuts[b]+numberOfRCuts[b]-1] = wolfCal.wolfCutoffCoulombEnd[b] * wolfCal.wolfCutoffCoulombEnd[b];
    }
    if (explicitlyAddEndRCut[b]){
      wolfAlpha[startOfNumAlphas[b]+numberOfAlphas[b]-1] = wolfCal.wolfAlphaEnd[b];
    }
    for(uint r = 1; r < numberOfRCuts[b]; r++) {
      for(uint a = 1 ; a < numberOfAlphas[b]; a++) {
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
    numberOfRCuts[b] += (int)((wolfCal.wolfCutoffCoulombEnd[b] - wolfCal.wolfCutoffCoulombStart[b]) / wolfCal.wolfCutoffCoulombDelta[b]);
    numberOfAlphas[b] += (int)((wolfCal.wolfAlphaEnd[b] - wolfCal.wolfAlphaStart[b]) / wolfCal.wolfAlphaDelta[b]);
    if (abs(wolfCal.wolfAlphaDelta[b] * numberOfAlphas[b] - wolfCal.wolfCutoffCoulombEnd[b]) > 0.01){
          numberOfAlphas[b] += 1;
          explicitlyAddEndAlpha[b] = true;
    } else {
          explicitlyAddEndAlpha[b] = false;
    }
    if (abs(wolfCal.wolfCutoffCoulombDelta[b] * numberOfRCuts[b] - wolfCal.wolfCutoffCoulombEnd[b]) > 1){
          numberOfRCuts[b] += 1;
          explicitlyAddEndRCut[b] = true;
    } else {
          explicitlyAddEndRCut[b] = false;
    }
  }
}

void WolfCalibration::AllocMem(){
  totalNumAlphas = 0;
  totalNumRCuts = 0;
  totalNumWolfFactors = 0;
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    startOfNumAlphas[b] += totalNumAlphas;
    startOfNumRCuts[b] += totalNumRCuts;
    startOfWolfFactor[b] += totalNumWolfFactors;
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