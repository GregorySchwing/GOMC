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
#include "System.h"
#include "CalculateEnergy.h"

class WolfCalibrationOutput : public OutputableBase
{
public:
  WolfCalibrationOutput(System & sys, StaticVals const& statV);

  ~WolfCalibrationOutput();

  virtual void DoOutput(const ulong step);
  virtual void DoOutputRestart(const ulong step) {}  
  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output);
  virtual void Sample(const ulong step) {}


private:

  void WriteHeader(void);
  void WriteGraceParFile(void);

  std::string GetString(double a, uint p);
  std::string GetString(ulong step);

  System & sysRef;
  StaticVals const& statValRef;
  const CalculateEnergy& calcEn;
  uint stepsPerSample;

  double ** electrostaticEnergies[BOX_TOTAL];

  //const CalculateEnergy& calcEn;
  std::ofstream outF[BOX_TOTAL];
  std::ofstream outFPar[BOX_TOTAL];
  std::string name[BOX_TOTAL];
  std::string namePar[BOX_TOTAL];
};

#endif