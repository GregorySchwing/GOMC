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
sysRef(sys), statValRef(statV)
{
      for(uint b = 0 ; b < BOX_TOTAL; b++) {
            electrostaticEnergies[b] =  new double*[statValRef.forcefield.numberOfRCuts[b]];
            for (int r = 0; r < statValRef.forcefield.numberOfRCuts[b]; ++r){
                  electrostaticEnergies[b][r] =  new double[statValRef.forcefield.numberOfAlphas[b]];
            }
      }
}

  WolfCalibrationOutput::~WolfCalibrationOutput()
  {
      for(uint b = 0 ; b < BOX_TOTAL; b++) {
            for (int r = 0; r < statValRef.forcefield.numberOfRCuts[b]; ++r){
                  delete[] electrostaticEnergies[b][r];
            }
            delete[] electrostaticEnergies[b];
      }
  }

void WolfCalibrationOutput::Init(pdb_setup::Atoms const& atoms,
                            config_setup::Output const& output) {
      stepsPerSample = output.wolfCalibration.settings.frequency;
      stepsPerOut = output.wolfCalibration.settings.frequency;
      enableOut = output.wolfCalibration.settings.enable;
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
            }
      //WriteHeader();
      }
}

void WolfCalibrationOutput::WriteHeader(void)
{
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    if (outF[b].is_open()) {
      std::string firstRow = "";
      std::string secondRow = "";
      firstRow += "RCutCoulomb ";
      for (int r = 0; r < statValRef.forcefield.numberOfRCuts[b]; ++r){
            firstRow += GetString(statValRef.forcefield.rCutCoulomb[b][r], 4);
            firstRow += ", ";

      }

      secondRow += "Alpha ";
      for (int a = 0; a < statValRef.forcefield.numberOfAlphas[b]; ++a){
            secondRow += GetString(statValRef.forcefield.wolfAlpha[b][a], 4);
            secondRow += ", ";
      }
      
      outF[b] << firstRow;
      outF[b] << std::endl;
      outF[b] << secondRow;
    } else {
      std::cerr << "Unable to write to file \"" <<  name[b] << "\" "
                << "(Wolf Calibration file)" << std::endl;
    }
  }
}

std::string WolfCalibrationOutput::GetString(double a, uint p)
{
  std::stringstream sstrm;
  std::string tempStr;
  sstrm << std::fixed << std::setprecision(p) << (a);
  //sstrm.precision(p);
  //sstrm >> tempStr;
  tempStr = sstrm.str();
  return tempStr;
}