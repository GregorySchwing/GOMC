
#include "EnsemblePreprocessor.h" //For BOX_TOTAL, etc.
#include "Setup.h"

class WolfCalibration{
public:
WolfCalibration(config_setup::WolfCalibration const& wolfCal);
~WolfCalibration();
double GetAlpha(int box, int indexForAlpha);
double GetRCut(int box, int indexForRCut);
double GetRCutSq(int box, int indexForRCut);
double GetWolfFactor1(int box, int indexForRCut, int indexForAlpha);
double GetWolfFactor2(int box, int indexForRCut, int indexForAlpha);
double GetWolfFactor3(int box, int indexForRCut, int indexForAlpha);
double GetNumberOfRCuts(int box);
double GetNumberOfAlphas(int box);
int GetStartOfWolfFactors(int box);
int GetTotalNumWolfFactors();

private:
  friend class WolfCalibrationOutput;
  void InitWolfCalibration(config_setup::WolfCalibration const& wolfCal);
  void CalculateWolfCalibrationMemoryUsage(config_setup::WolfCalibration const& wolfCal);
  void AllocMem();
  void DeallocMem();

  double * wolfAlpha; //alpha term for Wolf Electrostatic and constant factors
  double * wolfFactor1; //alpha term for Wolf Electrostatic and constant factors
  double * wolfFactor2;  //alpha term for Wolf Electrostatic and constant factors
  double * wolfFactor3; //alpha term for Wolf Electrostatic and constant factors
  double * rCutCoulomb;  //!<Cutoff Coulomb interaction(angstroms)
  double * rCutCoulombSq; //!<Cutoff Coulomb interaction(angstroms)
  int totalNumAlphas;
  int totalNumRCuts;
  int totalNumWolfFactors;
  int startOfNumRCuts[BOX_TOTAL];
  int startOfNumAlphas[BOX_TOTAL];
  int startOfWolfFactor[BOX_TOTAL];
  int numberOfRCuts[BOX_TOTAL];
  int numberOfAlphas[BOX_TOTAL];
  std::string wolfKindStrings[4] = {"HYBRID", "VLUGT", "GROSS", "VLUGTWINTRACUTOFF"};
  std::string coulKindStrings[2] = {"DSP", "DSF"};
  // Wolf Calibration
};