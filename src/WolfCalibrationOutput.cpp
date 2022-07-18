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
sysRef(sys), calcEn(sys.calcEnergy), statValRef(statV), wolfCalRef(statV.wolfCal)
{
      // This is neccessary to check for correctness of single point energy calculations.
      printOnFirstStep = true;
}

WolfCalibrationOutput::~WolfCalibrationOutput()
{

}

void WolfCalibrationOutput::Init(pdb_setup::Atoms const& atoms,
                            config_setup::Output const& output) {
      for (int b = 0; b < BOX_TOTAL; ++b){
            numberOfRCuts[b] = wolfCalRef.GetNumberOfRCuts(b);
            numberOfAlphas[b] = wolfCalRef.GetNumberOfAlphas(b);
            startOfWolfFactor[b] = wolfCalRef.GetStartOfWolfFactors(b);
      }
      stepsPerSample = output.wolfCalibration.settings.frequency;
      stepsPerOut = output.wolfCalibration.settings.frequency;
      enableOut = output.wolfCalibration.settings.enable;
      electrostaticEnergies.resize(BOX_TOTAL*WOLF_TOTAL_KINDS*COUL_TOTAL_KINDS*wolfCalRef.GetTotalNumWolfFactors());
      electrostaticEnergies.assign(BOX_TOTAL*WOLF_TOTAL_KINDS*COUL_TOTAL_KINDS*wolfCalRef.GetTotalNumWolfFactors(), 0.0);
      if(enableOut) {
            for (uint b = 0; b < BOX_TOTAL; ++b) {
                  for (uint wolfKind = 0; wolfKind < WOLF_TOTAL_KINDS; ++wolfKind){
                        for (uint coulKind = 0; coulKind < COUL_TOTAL_KINDS; ++coulKind){
                              std::stringstream sstrm;
                              std::string strKind, fileName;
                              sstrm << (b);
                              sstrm >> strKind;
                              fileName = "Wolf_Calibration_";
                              fileName += wolfCalRef.wolfKindStrings[wolfKind];
                              fileName += "_";
                              fileName += wolfCalRef.coulKindStrings[coulKind];
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
            for (int r = 0; r < wolfCalRef.numberOfRCuts[b]; ++r){
                  for (int a = 0; a < wolfCalRef.numberOfAlphas[b]; ++a){
                        firstRow += "(";
                        firstRow += GetString(wolfCalRef.GetRCut(b, r), 4);
                        firstRow += ", ";
                        firstRow += GetString(wolfCalRef.GetAlpha(b, a), 4);
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
            for (int r = 0; r < wolfCalRef.numberOfRCuts[b]; ++r){
                  for (int a = 0; a < wolfCalRef.numberOfAlphas[b]; ++a){
                        firstRow += "\ts";
                        firstRow += GetString(counter);
                        firstRow += " legend \"(";
                        firstRow += GetString(wolfCalRef.rCutCoulomb[wolfCalRef.startOfNumRCuts[b]+r], 4);
                        firstRow += ", ";
                        firstRow += GetString(wolfCalRef.wolfAlpha[wolfCalRef.startOfNumAlphas[b]+a], 4);
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
      calcEn.WolfCalibrationEnergy(&electrostaticEnergies[0]);

      // Calc reference epot
      sysRef.SwapWolfAndEwaldPointers();
      
      statValRef.forcefield.ewald = true;
      statValRef.forcefield.wolf = false;
      // Set isVlugtWolf == false
      // Since ewald and wolf share a ff object
      statValRef.forcefield.SetWolfKind(0);
      SystemPotential ewaldRef = calcEn.SystemTotal();
      ewaldRef.Total();
      printf("ew en %f\n", ewaldRef.boxEnergy[0].total);
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
                        for (int r = 0; r < wolfCalRef.numberOfRCuts[box]; ++r){
                              for (int a = 0; a < wolfCalRef.numberOfAlphas[box]; ++a){
                                    // If you dont use std::abs, double is converted to int 
                                    row += GetString((std::abs(ewaldRef.boxEnergy[box].total) -  std::abs(electrostaticEnergies[wolfCalRef.GetIndex(box, wolfKind, coulKind, r, a)]))/ std::abs(ewaldRef.boxEnergy[box].total), 8);
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