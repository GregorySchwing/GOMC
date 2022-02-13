/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#ifndef WOLF_CALIBRATION_OUTPUT_H
#define WOLF_CALIBRATION_OUTPUT_H

#include "OutputAbstracts.h"
#include <iostream>
#include "GOMC_Config.h"


class WolfCalibrationOutput : public OutputableBase
{
public:
  WolfCalibrationOutput(System & sys, StaticVals const& statV);

  ~WolfCalibrationOutput();

  virtual void DoOutput(const ulong step) {}
  virtual void DoOutputRestart(const ulong step) {}  
  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output);
  virtual void Sample(const ulong step) {}


private:

  //void WriteHeader(void);
  std::string GetString(double a, uint p);

  System & sysRef;
  StaticVals const& statValRef;
  uint stepsPerSample;
  const config_setup::WolfCalibration&  wolfCal;

  uint numberOfRCutStates[BOX_TOTAL];
  uint numberOfAlphaStates[BOX_TOTAL];

  double * alphas;
  double * rcutcoulombs;
  double * electrostaticEnergies;


  bool explicitlyAddEndAlpha[BOX_TOTAL];
  bool explicitlyAddEndRCut[BOX_TOTAL];

  //const CalculateEnergy& calcEn;
  std::ofstream outF[BOX_TOTAL];
  std::string name[BOX_TOTAL];
};

#endif