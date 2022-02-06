/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include <stdint.h>
#include "WolfCalibrationOutput.h"
#include "GOMC_Config.h"


WolfCalibrationOutput::WolfCalibrationOutput(System & sys, StaticVals const& statV):
sysRef(sys), statValRef(statV), wolfCal(statV.wolfCal)
{
}

void WolfCalibrationOutput::Init(pdb_setup::Atoms const& atoms,
                            config_setup::Output const& output) {
      stepsPerSample = output.wolfCalibration.settings.frequency;
      stepsPerOut = output.wolfCalibration.settings.frequency;
      enableOut = output.wolfCalibration.settings.enable;
      // Still need to implement these.
      for (uint b = 0; b < BOX_TOTAL; ++b) {
            numberOfRCutStates[b] = (int)((wolfCal.wolfCutoffCoulombEnd[b] - wolfCal.wolfCutoffCoulombStart[b]) / wolfCal.wolfCutoffCoulombDelta[b]);
            numberOfAlphaStates[b] = (int)((wolfCal.wolfAlphaEnd[b] - wolfCal.wolfAlphaStart[b]) / wolfCal.wolfAlphaDelta[b]);
      }
      if(enableOut) {
            for (uint b = 0; b < BOX_TOTAL; ++b) {
                  std::stringstream sstrm;
                  std::string strKind, fileName;
                  sstrm << (b);
                  sstrm >> strKind;
                  fileName = "Wolf_Calibration_BOX_";
                  fileName += strKind;
                  fileName += "_";
                  fileName += uniqueName;
                  fileName += ".dat";
                  #if GOMC_LIB_MPI
                        name[b] = pathToReplicaOutputDirectory + fileName;
                  #else
                        name[b] = fileName;
                  #endif
                  outF[b].open(name[b].c_str(), std::ofstream::out);
                  //energyDiff[b] = new Energy[lambdaSize];
            }
      //WriteHeader();
      }
}