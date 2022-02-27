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
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                        electrostaticEnergies[b][wolfKind][coulKind] =  new double*[statValRef.forcefield.numberOfRCuts[b]];
                        for (int r = 0; r < statValRef.forcefield.numberOfRCuts[b]; ++r){
                              electrostaticEnergies[b][wolfKind][coulKind][r] =  new double[statValRef.forcefield.numberOfAlphas[b]];
                        }
                  }
            }
      }
}

  WolfCalibrationOutput::~WolfCalibrationOutput()
  {
      for(uint b = 0 ; b < BOX_TOTAL; b++) {
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                        for (int r = 0; r < statValRef.forcefield.numberOfRCuts[b]; ++r){
                              delete[] electrostaticEnergies[b][wolfKind][coulKind][r];
                        }
                        delete[] electrostaticEnergies[b][wolfKind][coulKind];
                  }
            }
      }
  }

void WolfCalibrationOutput::Init(pdb_setup::Atoms const& atoms,
                            config_setup::Output const& output) {
      stepsPerSample = output.wolfCalibration.settings.frequency;
      stepsPerOut = output.wolfCalibration.settings.frequency;
      enableOut = output.wolfCalibration.settings.enable;
      if(enableOut) {
            for (uint b = 0; b < BOX_TOTAL; ++b) {
                  for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                        for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                              std::stringstream sstrm;
                              std::string strKind, fileName;
                              sstrm << (b);
                              sstrm >> strKind;
                              fileName = "Wolf_Calibration_";
                              fileName += statValRef.forcefield.wolfKindStrings[wolfKind];
                              fileName += "_";
                              fileName += statValRef.forcefield.coulKindStrings[coulKind];
                              fileName += "_BOX_";
                              fileName += strKind;
                              fileName += "_";
                              fileName += uniqueName;
                              #if GOMC_LIB_MPI
                                    name[b][wolfKind][coulKind] = pathToReplicaOutputDirectory + fileName + ".dat";
                                    namePar[b][wolfKind][coulKind] = pathToReplicaOutputDirectory + fileName + ".par";
                              #else
                                    name[b][wolfKind][coulKind] = fileName + ".dat";
                                    namePar[b][wolfKind][coulKind] = fileName + ".par";

                              #endif
                              outF[b][wolfKind][coulKind].open(name[b][wolfKind][coulKind].c_str(), std::ofstream::out);
                              outFPar[b][wolfKind][coulKind].open(namePar[b][wolfKind][coulKind].c_str(), std::ofstream::out);
                              WriteHeader(b, wolfKind, coulKind);
                              WriteGraceParFile(b, wolfKind, coulKind);
                        }
                  }
            }
      }
}

void WolfCalibrationOutput::WriteHeader(uint b, uint wolfKind, uint coulKind)
{
      if (outF[b][wolfKind][coulKind].is_open()) {
            std::string firstRow = "";
            firstRow += "Step#\t";
            // We skip the reference r cut with reference alpha.
            // r = 0, a = 0
            // So there are no duplicate columns.
            for (int r = 1; r < statValRef.forcefield.numberOfRCuts[b]; ++r){
                  for (int a = 1; a < statValRef.forcefield.numberOfAlphas[b]; ++a){
                        firstRow += "(";
                        firstRow += GetString(statValRef.forcefield.rCutCoulomb[b][r], 4);
                        firstRow += ", ";
                        firstRow += GetString(statValRef.forcefield.wolfAlpha[b][a], 4);
                        firstRow += ")\t";
                  }
            }
            outF[b][wolfKind][coulKind] << firstRow;
            outF[b][wolfKind][coulKind] << std::endl;
      } else {
            std::cerr << "Unable to write to file \"" <<  name[b][wolfKind][coulKind] << "\" "
                        << "(Wolf Calibration file)" << std::endl;
      }
}


void WolfCalibrationOutput::WriteGraceParFile(uint b, uint wolfKind, uint coulKind)
{
      for (uint b = 0; b < BOX_TOTAL; ++b) {
            if (outFPar[b][wolfKind][coulKind].is_open()) {
                  int counter = 0;
                  std::string firstRow = "";
                  firstRow += "with g0\n";
                  // We skip the reference r cut with reference alpha.
                  // r = 0, a = 0
                  // So there are no duplicate columns.
                  for (int r = 1; r < statValRef.forcefield.numberOfRCuts[b]; ++r){
                        for (int a = 1; a < statValRef.forcefield.numberOfAlphas[b]; ++a){
                              firstRow += "\ts";
                              firstRow += GetString(counter);
                              firstRow += " legend \"(";
                              firstRow += GetString(statValRef.forcefield.rCutCoulomb[b][r], 4);
                              firstRow += ", ";
                              firstRow += GetString(statValRef.forcefield.wolfAlpha[b][a], 4);
                              firstRow += ")\"\n";
                              ++counter;
                        }
                  }
                  outFPar[b][wolfKind][coulKind] << firstRow;
                  outFPar[b][wolfKind][coulKind] << std::endl;
            } else {
                  std::cerr << "Unable to write to file \"" <<  name[b] << "\" "
                              << "(Wolf Calibration file)" << std::endl;
            }
      }
}

void WolfCalibrationOutput::DoOutput(const ulong step) {
      uint wolfKindOrig = sysRef.calcEwald->GetWolfKind();
      uint coulKindOrig = sysRef.calcEwald->GetCoulKind();
      calcEn.WolfCalibrationEnergy(electrostaticEnergies);
      sysRef.calcEwald->SetWolfKind(wolfKindOrig);
      sysRef.calcEwald->SetCoulKind(coulKindOrig);
      // Eventually use this to calc refernce
      statValRef.forcefield.ewald = true;
      sysRef.SwapWolfAndEwaldPointers();
      SystemPotential ewaldRef = calcEn.SystemTotal();
      statValRef.forcefield.ewald = false;
      sysRef.SwapWolfAndEwaldPointers();
      std::string row = "";
      row += GetString(step);
      row += "\t";

      for (uint box = 0; box < BOX_TOTAL; ++box) {       
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){           
                        // We skip the reference r cut with reference alpha.
                        // r = 0, a = 0
                        // So there are no duplicate columns.
                        for (int r = 1; r < statValRef.forcefield.numberOfRCuts[box]; ++r){
                              for (int a = 1; a < statValRef.forcefield.numberOfAlphas[box]; ++a){
                                    row += GetString((abs(ewaldRef.boxEnergy[box].total) -  abs(electrostaticEnergies[box][wolfKind][coulKind][r][a]))/ abs(ewaldRef.boxEnergy[box].total), 4);
                                    row += "\t";
                              }
                        }
                        outF[box][wolfKind][coulKind] << row;
                        outF[box][wolfKind][coulKind] << std::endl;
                  }
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

std::string WolfCalibrationOutput::GetString(ulong step)
{
      std::stringstream ss;
      ss << step;
      return ss.str();
}