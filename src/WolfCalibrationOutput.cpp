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

  WolfCalibrationOutput::~WolfCalibrationOutput()
  {
      delete[] alphas;
      delete[] rcutcoulombs;
      delete[] energyDiff;
  }

void WolfCalibrationOutput::Init(pdb_setup::Atoms const& atoms,
                            config_setup::Output const& output) {
      stepsPerSample = output.wolfCalibration.settings.frequency;
      stepsPerOut = output.wolfCalibration.settings.frequency;
      enableOut = output.wolfCalibration.settings.enable;
      int totalRCutStates = 0;
      int totalAlphaStates= 0;
      int totaEnergyDiffStates = 0;
      // Calculate the number of calibration points from the ranges provided
      // If delta is not a common multiple of the (End - Start) explicitly add the End.
      for (uint b = 0; b < BOX_TOTAL; ++b) {
            numberOfRCutStates[b] = (int)((wolfCal.wolfCutoffCoulombEnd[b] - wolfCal.wolfCutoffCoulombStart[b]) / wolfCal.wolfCutoffCoulombDelta[b]);
            numberOfAlphaStates[b] = (int)((wolfCal.wolfAlphaEnd[b] - wolfCal.wolfAlphaStart[b]) / wolfCal.wolfAlphaDelta[b]);
            if (abs(wolfCal.wolfAlphaDelta[b] * numberOfAlphaStates[b] - wolfCal.wolfCutoffCoulombEnd[b]) > 0.01){
                  numberOfAlphaStates[b]= numberOfAlphaStates[b] + 1;
                  explicitlyAddEndAlpha[b] = true;
            } else {
                  explicitlyAddEndAlpha[b] = false;
            }
            if (abs(wolfCal.wolfCutoffCoulombDelta[b] * numberOfRCutStates[b] - wolfCal.wolfCutoffCoulombEnd[b]) > 1){
                  numberOfRCutStates[b]= numberOfRCutStates[b] + 1;
                  explicitlyAddEndRCut[b] = true;
            } else {
                  explicitlyAddEndRCut[b] = false;
            }
            totalRCutStates += numberOfRCutStates[b];
            totalAlphaStates += numberOfAlphaStates[b];
            totaEnergyDiffStates += numberOfRCutStates[b] * numberOfAlphaStates[b];
      }
      rcutcoulombs = new double[totalRCutStates];
      alphas = new double[totalAlphaStates];
      energyDiff =  new Energy[totaEnergyDiffStates];
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
/*
void WolfCalibrationOutput::WriteHeader(void)
{
  for (uint b = 0; b < BOX_TOTAL; ++b) {
    if (outF[b].is_open()) {
      std::string firstRow = "";
      std::string secondRow = "";
      toPrint += "Alpha ";
      for (int i = 0; i < numberOfRCutStates[b]; ++i){
            for (int i = 0; i < numberOfAlphaStates[b]; ++i){

            }
      }
      toPrint += GetString(var->T_in_K, 4);
      toPrint += "(K), Lambda State ";
      toPrint += GetString(iState, 0);
      toPrint += ": (lambda Coulomb, lambda VDW) = (";
      toPrint += GetString(freeEnVal.lambdaCoulomb[iState], 4);
      toPrint += ",";
      toPrint += GetString(freeEnVal.lambdaVDW[iState], 4);
      toPrint += ")\n";
      outF[b] << toPrint;

      //We care about format
      outF[b] << std::setw(11) << std::left << "#Steps" << " ";
      outF[b] << std::setw(25) << std::right << "Total_En(kJ/mol)" << " ";
      toPrint = "dU/dL(Coulomb=";
      toPrint += GetString(freeEnVal.lambdaCoulomb[iState], 4);
      toPrint += ")";
      outF[b] << std::setw(25) << std::right << toPrint << " ";
      toPrint = "dU/dL(VDW=";
      toPrint += GetString(freeEnVal.lambdaVDW[iState], 4);
      toPrint += ")";
      outF[b] << std::setw(25) << std::right << toPrint << " ";

      std::string fixStr = "DelE(L->(";
      for(uint i = 0; i < lambdaSize; i++) {
        toPrint = fixStr;
        toPrint += GetString(freeEnVal.lambdaCoulomb[i], 4);
        toPrint += ",";
        toPrint += GetString(freeEnVal.lambdaVDW[i], 4);
        toPrint += "))";
        outF[b] << std::setw(25) << std::right << toPrint << " ";
      }
#if ENSEMBLE == NVT
      if(var->pressureCalc) {
        outF[b] << std::setw(25) << std::right << "PV(kJ/mol)";
      }
#elif ENSEMBLE == NPT
      outF[b] << std::setw(25) << std::right << "PV(kJ/mol)";
#endif
      outF[b] << std::endl;
      outF[b] << std::setprecision(10);
      outF[b].setf(std::ios_base::right, std::ios_base::adjustfield);
    } else
      std::cerr << "Unable to write to file \"" <<  name[b] << "\" "
                << "(Free Energy file)" << std::endl;
  }
}
*/
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