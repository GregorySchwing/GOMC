/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include <stdint.h>
#include "WolfCalibrationOutput.h"
#include "GOMC_Config.h"


WolfCalibrationOutput::WolfCalibrationOutput(System & sys, StaticVals & statV):
sysRef(sys), calcEn(sys.calcEnergy), statValRef(statV)
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
                  #if GOMC_LIB_MPI
                        name[b] = pathToReplicaOutputDirectory + fileName + ".dat";
                        namePar[b] = pathToReplicaOutputDirectory + fileName + ".par";
                  #else
                        name[b] = fileName + ".dat";
                        namePar[b] = fileName + ".par";

                  #endif
                  outF[b].open(name[b].c_str(), std::ofstream::out);
                  outFPar[b].open(namePar[b].c_str(), std::ofstream::out);
            }
            WriteHeader();
            WriteGraceParFile();
      }
}

void WolfCalibrationOutput::WriteHeader(void)
{
      for (uint b = 0; b < BOX_TOTAL; ++b) {
            if (outF[b].is_open()) {
                  std::string firstRow = "";
                  firstRow += "Step#\t";
                  for (int r = 0; r < statValRef.forcefield.numberOfRCuts[b]; ++r){
                        for (int a = 0; a < statValRef.forcefield.numberOfAlphas[b]; ++a){
                              firstRow += "(";
                              firstRow += GetString(statValRef.forcefield.rCutCoulomb[b][r], 4);
                              firstRow += ", ";
                              firstRow += GetString(statValRef.forcefield.wolfAlpha[b][a], 4);
                              firstRow += ")\t";
                              // We only want the reference r cut with reference alpha.
                              if (r == 0)
                                    break;
                        }
                  }
                  outF[b] << firstRow;
                  outF[b] << std::endl;
            } else {
                  std::cerr << "Unable to write to file \"" <<  name[b] << "\" "
                              << "(Wolf Calibration file)" << std::endl;
            }
      }
}


void WolfCalibrationOutput::WriteGraceParFile(void)
{
      for (uint b = 0; b < BOX_TOTAL; ++b) {
            if (outFPar[b].is_open()) {
                  int counter = 0;
                  std::string firstRow = "";
                  firstRow += "with g0\n";
                  for (int r = 0; r < statValRef.forcefield.numberOfRCuts[b]; ++r){
                        for (int a = 0; a < statValRef.forcefield.numberOfAlphas[b]; ++a){
                              firstRow += "\ts";
                              firstRow += GetString(counter);
                              firstRow += " legend \"(";
                              firstRow += GetString(statValRef.forcefield.rCutCoulomb[b][r], 4);
                              firstRow += ", ";
                              firstRow += GetString(statValRef.forcefield.wolfAlpha[b][a], 4);
                              firstRow += ")\"\n";
                              ++counter;
                              // We only want the reference r cut with reference alpha.
                              if (r == 0)
                                    break;
                        }
                  }
                  outFPar[b] << firstRow;
                  outFPar[b] << std::endl;
            } else {
                  std::cerr << "Unable to write to file \"" <<  name[b] << "\" "
                              << "(Wolf Calibration file)" << std::endl;
            }
      }
}

void WolfCalibrationOutput::DoOutput(const ulong step) {
      for (uint box = 0; box < BOX_TOTAL; ++box) {
            calcEn.WolfCalibrationEnergyChange(box,
                                          electrostaticEnergies);
            statValRef.forcefield.ewald = true;
            std::string row = "";
            row += GetString(step);
            row += "\t";
            for (int r = 0; r < statValRef.forcefield.numberOfRCuts[box]; ++r){
                  for (int a = 0; a < statValRef.forcefield.numberOfAlphas[box]; ++a){
                        row += GetString(electrostaticEnergies[box][r][a], 4);
                        row += "\t";
                        // We only want the reference r cut with reference alpha.
                        if (r == 0)
                              break;
                  }
            }
            outF[box] << row;
            outF[box] << std::endl;
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

std::string WolfCalibrationOutput::GetString(ulong step)
{
      std::stringstream ss;
      ss << step;
      return ss.str();
}