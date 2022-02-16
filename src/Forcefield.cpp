/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Forcefield.h" //Header spec.
//Setup partner classes
#include "Setup.h"
#include "FFShift.h"
#include "FFSwitch.h"
#include "FFSwitchMartini.h"
#include "FFExp6.h"
#define _USE_MATH_DEFINES
#include <cmath>

Forcefield::Forcefield()
{
  particles = NULL;
  angles = NULL;
  OneThree = false; //default behavior is to turn off 1-3 interaction
  OneFour = true;   // to turn on 1-4 interaction
  OneN = true;      // and turn on 1-n interaction
  // Default, when not calibrating wolf, these values are 1.
  for (int box = 0; box < BOX_TOTAL; ++box){
    numberOfRCuts[box] = 1;
    numberOfAlphas[box] = 1;
  }
}

Forcefield::~Forcefield()
{
  if(particles != NULL)
    delete particles;
  if( angles != NULL)
    delete angles;
  DeallocMem();
}

void Forcefield::Init(const Setup& set,
                      config_setup::WolfCalibration const& wolfCal)
{
  wolfCalibration = wolfCal.enable;
  if(wolfCalibration){
    CalculateWolfCalibrationMemoryUsage(wolfCal);
  }
  AllocMem();
  InitBasicVals(set.config.sys, set.config.in.ffKind);
  if(wolfCalibration){
    InitWolfCalibration(wolfCal);
  }
  particles->Init(set.ff.mie, set.ff.nbfix);
  bonds.Init(set.ff.bond);
  angles->Init(set.ff.angle);
  dihedrals.Init(set.ff.dih);
  coulKind = set.config.sys.ff.COUL_KIND;
  wolfKind = set.config.sys.ff.WOLF_KIND;
  // Only Vlugt Wolf alters the FF behavior by removing cutoffs for Intra Undampened
  if (wolfKind == 1){
    isVlugtWolf = true;
  } else {
    isVlugtWolf = false;
  }
}

void Forcefield::AllocMem(){
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    wolfAlpha[b] = new double[numberOfAlphas[b]];
    rCutCoulomb[b] = new double[numberOfRCuts[b]];
    rCutCoulombSq[b] = new double[numberOfRCuts[b]];
    wolfFactor1[b] = new double*[numberOfRCuts[b]];
    wolfFactor2[b] = new double*[numberOfRCuts[b]];
    wolfFactor3[b] = new double*[numberOfRCuts[b]];
    for (int r = 0; r < numberOfRCuts[b]; ++r){
      wolfFactor1[b][r] = new double[numberOfAlphas[b]];
      wolfFactor2[b][r] = new double[numberOfAlphas[b]];
      wolfFactor3[b][r] = new double[numberOfAlphas[b]];
    }
  }
}

void Forcefield::DeallocMem(){
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    delete[] wolfAlpha[b];
    delete[] rCutCoulomb[b];
    delete[] rCutCoulombSq[b];
    for (int r = 0; r < numberOfRCuts[b]; ++r){
      delete[] wolfFactor1[b][r];
      delete[] wolfFactor2[b][r];
      delete[] wolfFactor3[b][r];
    }
    delete[] wolfFactor1[b];
    delete[] wolfFactor2[b];
    delete[] wolfFactor3[b];
  }
}

void Forcefield::InitBasicVals(config_setup::SystemVals const& val,
                               config_setup::FFKind const& ffKind)
{
  useLRC = val.ff.doTailCorr;
  T_in_K = val.T.inKelvin;
  rCut = val.ff.cutoff;
  rCutSq = rCut * rCut;
  rCutLow = val.ff.cutoffLow;
  rCutLowSq = rCutLow * rCutLow;
  scaling_14 = val.elect.oneFourScale;
  beta = 1 / T_in_K;

  vdwKind = val.ff.VDW_KIND;
  exckind = val.exclude.EXCLUDE_KIND;
  freeEnergy = val.freeEn.enable;

  electrostatic = val.elect.enable;
  ewald = val.elect.ewald;
  wolf = val.elect.wolf;
  tolerance = val.elect.tolerance;
  rswitch = val.ff.rswitch;
  dielectric = val.elect.dielectric;

  if(val.freeEn.enable) {
    sc_alpha = val.freeEn.scaleAlpha;
    sc_sigma = val.freeEn.scaleSigma;
    sc_power = val.freeEn.scalePower;
    sc_coul = val.freeEn.scaleCoulomb;
  } else if (val.neMTMCVal.enable) {
    sc_alpha = val.neMTMCVal.scaleAlpha;
    sc_sigma = val.neMTMCVal.scaleSigma;
    sc_power = val.neMTMCVal.scalePower;
    sc_coul = val.neMTMCVal.scaleCoulomb;
  } else {
    sc_alpha = 0.0;
    sc_sigma = 0.0;
    sc_power = 0;
    sc_coul = false;
  }
  sc_sigma_6 = pow(sc_sigma, 6.0);

  for(uint b = 0 ; b < BOX_TOTAL; b++) {
    rCutCoulomb[b][0] = val.elect.cutoffCoulomb[b];
    rCutCoulombSq[b][0] = rCutCoulomb[b][0] * rCutCoulomb[b][0];
    alpha[b] = sqrt(-log(tolerance)) / rCutCoulomb[b][0];
    alphaSq[b] = alpha[b] * alpha[b];
    recip_rcut[b] = -2.0 * log(tolerance) / rCutCoulomb[b][0];
    recip_rcut_Sq[b] = recip_rcut[b] * recip_rcut[b];
    if (wolf){
      wolfAlpha[b][0] = val.elect.wolfAlpha[b];
      wolfFactor1[b][0][0] = erfc(wolfAlpha[b][0]*rCutCoulomb[b][0])/rCutCoulomb[b][0];
      wolfFactor2[b][0][0] = wolfFactor1[b][0][0]/rCutCoulomb[b][0];
      wolfFactor2[b][0][0] += wolfAlpha[b][0] *  M_2_SQRTPI * 
                        exp(-1.0*wolfAlpha[b][0]*wolfAlpha[b][0]*rCutCoulombSq[b][0])
                        /rCutCoulomb[b][0];
      wolfFactor3[b][0][0] = wolfAlpha[b][0] *  M_2_SQRTPI;
    }
  }

  vdwGeometricSigma = val.ff.vdwGeometricSigma;
  isMartini = ffKind.isMARTINI;
  exp6 = (vdwKind == val.ff.VDW_EXP6_KIND);

#if ENSEMBLE == GCMC
  isFugacity = val.chemPot.isFugacity;
#endif

  if(vdwKind == val.ff.VDW_STD_KIND)
    particles = new FFParticle(*this);
  else if(vdwKind == val.ff.VDW_EXP6_KIND)
    particles = new FF_EXP6(*this);
  else if(vdwKind == val.ff.VDW_SHIFT_KIND)
    particles = new FF_SHIFT(*this);
  else if (vdwKind == val.ff.VDW_SWITCH_KIND && ffKind.isMARTINI)
    particles = new FF_SWITCH_MARTINI(*this);
  else if (vdwKind == val.ff.VDW_SWITCH_KIND && !ffKind.isMARTINI)
    particles = new FF_SWITCH(*this);
  else {
    std::cout << "Undefined Potential Type detected!\n" << "Exiting!\n";
    exit(EXIT_FAILURE);
  }


  if(ffKind.isMARTINI)
    angles = new FFAngleMartini();
  else
    angles = new FFAngles();

  // Define type of interaction to be included. ex. 1-3, 1-4 and more
  if(exckind == val.exclude.EXC_ONETWO_KIND) {
    OneThree = true, OneFour = true, OneN = true;
  } else if(exckind == val.exclude.EXC_ONETHREE_KIND) {
    OneThree = false, OneFour = true, OneN = true;
  } else if(exckind == val.exclude.EXC_ONEFOUR_KIND) {
    OneThree = false, OneFour = false, OneN = true;
  } else if(exckind == val.exclude.EXC_ONEN_KIND) {
    OneThree = false, OneFour = false, OneN = false;
  } else {
    std::cout << "Error: Unknown exclude value.\n";
    exit(EXIT_FAILURE);
  }

}

void Forcefield::CalculateWolfCalibrationMemoryUsage(config_setup::WolfCalibration const& wolfCal){
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

void Forcefield::InitWolfCalibration(config_setup::WolfCalibration const& wolfCal){
  for(uint b = 0 ; b < BOX_TOTAL; b++) {
    // Start at 1, since 0th index is from the config file and initted in InitBasicVals
    for(uint r = 0; r < numberOfRCuts[b]; r++) {
      // Start at 1, since 0th index is from the config file and initted in InitBasicVals
      rCutCoulomb[b][r+1] = wolfCal.wolfCutoffCoulombStart[b] + r*wolfCal.wolfCutoffCoulombDelta[b];
      rCutCoulombSq[b][r+1] = rCutCoulomb[b][r+1] * rCutCoulomb[b][r+1];
    }
    for(uint a = 0 ; a < numberOfAlphas[b]; a++) {
      wolfAlpha[b][a+1] = wolfCal.wolfAlphaStart[b] + a*wolfCal.wolfAlphaDelta[b];
    }
    for(uint r = 1; r < numberOfRCuts[b]; r++) {
      for(uint a = 1 ; a < numberOfAlphas[b]; a++) {
        wolfFactor1[b][r][a] = erfc(wolfAlpha[b][a]*rCutCoulomb[b][r])/rCutCoulomb[b][r];
        wolfFactor2[b][r][a] = wolfFactor1[b][r][a]/rCutCoulomb[b][r];
        wolfFactor2[b][r][a] += wolfAlpha[b][a] *  M_2_SQRTPI * 
                          exp(-1.0*wolfAlpha[b][a]*wolfAlpha[b][a]*rCutCoulombSq[b][r])
                          /rCutCoulomb[b][r];
        wolfFactor3[b][r][a] = wolfAlpha[b][a] *  M_2_SQRTPI;
      }
    }
  }
}

