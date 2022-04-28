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
sysRef(sys), calcEn(sys.calcEnergy), statValRef(statV), wolfCal(statV.forcefield.wolfCal)
{

}

  WolfCalibrationOutput::~WolfCalibrationOutput()
  {

  }

void WolfCalibrationOutput::Init(pdb_setup::Atoms const& atoms,
                            config_setup::Output const& output) {
      stepsPerSample = output.wolfCalibration.settings.frequency;
      stepsPerOut = output.wolfCalibration.settings.frequency;
      enableOut = output.wolfCalibration.settings.enable;
      for(uint b = 0 ; b < BOX_TOTAL; b++) {
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                        electrostaticEnergies[b][wolfKind][coulKind] =  new double*[wolfCal->numberOfRCuts[b]];
                        for (int r = 0; r < wolfCal->numberOfRCuts[b]; ++r){
                              electrostaticEnergies[b][wolfKind][coulKind][r] =  new double[wolfCal->numberOfAlphas[b]];
                        }
                  }
            }
      }
      if(enableOut) {
            for (uint b = 0; b < BOX_TOTAL; ++b) {
                  for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                        for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                              std::stringstream sstrm;
                              std::string strKind, fileName;
                              sstrm << (b);
                              sstrm >> strKind;
                              fileName = "Wolf_Calibration_";
                              fileName += wolfCal->wolfKindStrings[wolfKind];
                              fileName += "_";
                              fileName += wolfCal->coulKindStrings[coulKind];
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
            for (int r = 1; r < wolfCal->numberOfRCuts[b]; ++r){
                  for (int a = 1; a < wolfCal->numberOfAlphas[b]; ++a){
                        firstRow += "(";
                        firstRow += GetString(wolfCal->rCutCoulomb[wolfCal->startOfNumRCuts[b]+r], 4);
                        firstRow += ", ";
                        firstRow += GetString(wolfCal->wolfAlpha[wolfCal->startOfNumAlphas[b]+a], 4);
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
      if (outFPar[b][wolfKind][coulKind].is_open()) {
            int counter = 0;
            std::string firstRow = "";
            firstRow += "with g0\n";
            // We skip the reference r cut with reference alpha.
            // r = 0, a = 0
            // So there are no duplicate columns.
            for (int r = 1; r < wolfCal->numberOfRCuts[b]; ++r){
                  for (int a = 1; a < wolfCal->numberOfAlphas[b]; ++a){
                        firstRow += "\ts";
                        firstRow += GetString(counter);
                        firstRow += " legend \"(";
                        firstRow += GetString(wolfCal->rCutCoulomb[wolfCal->startOfNumRCuts[b]+r], 4);
                        firstRow += ", ";
                        firstRow += GetString(wolfCal->wolfAlpha[wolfCal->startOfNumAlphas[b]+a], 4);
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

void WolfCalibrationOutput::DoOutput(const ulong step) {
      uint wolfKindOrig = sysRef.calcEwald->GetWolfKind();
      uint coulKindOrig = sysRef.calcEwald->GetCoulKind();
      //calcEn.WolfCalibrationEnergy(electrostaticEnergies);

      // Calc reference epot
      sysRef.SwapWolfAndEwaldPointers();
      
      statValRef.forcefield.ewald = true;
      statValRef.forcefield.wolf = false;
      // Set isVlugtWolf == false
      // Since ewald and wolf share a ff object
      statValRef.forcefield.SetWolfKind(0);
      SystemPotential ewaldRef = calcEn.SystemTotal();

      sysRef.SwapWolfAndEwaldPointers();

      // Restore original wolf settings
      statValRef.forcefield.ewald = false;
      statValRef.forcefield.wolf = true;
      
      // Restore inter wolf settings
      statValRef.forcefield.SetWolfKind(wolfKindOrig);
      statValRef.forcefield.SetCoulKind(coulKindOrig);

      // Restore intra wolf settings
      sysRef.calcEwald->SetWolfKind(wolfKindOrig);
      sysRef.calcEwald->SetCoulKind(coulKindOrig);
      std::string row = "";
      row += GetString(step);
      row += "\t";

      for (uint box = 0; box < BOX_TOTAL; ++box) {       
            for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                  for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){  
                        std::string row = "";
                        row += GetString(step);
                        row += "\t";
                        // We skip the reference r cut with reference alpha.
                        // r = 0, a = 0
                        // So there are no duplicate columns.
                        for (int r = 1; r < wolfCal->numberOfRCuts[box]; ++r){
                              for (int a = 1; a < wolfCal->numberOfAlphas[box]; ++a){
                                    // If you dont use std::abs, double is converted to int 
                                    row += GetString((std::abs(ewaldRef.boxEnergy[box].total) -  std::abs(electrostaticEnergies[box][wolfKind][coulKind][r][a]))/ std::abs(ewaldRef.boxEnergy[box].total), 8);
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