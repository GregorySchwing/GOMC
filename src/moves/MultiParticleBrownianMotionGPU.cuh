/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.50
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in License.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef MULTIPARTICLEBROWNIANMOTIONGPU_H
#define MULTIPARTICLEBROWNIANMOTIONGPU_H

#include "MultiParticleBrownianMotion.h"

class MultiParticleBrownianGPU : public MultiParticleBrownian
{
public:
  MultiParticleBrownianGPU(System &sys, StaticVals const& statV);
  ~MultiParticleBrownianGPU() {

  }
/*
  virtual uint Prep(const double subDraw, const double movPerc);
  // To relax the system in NE_MTMC move
  virtual uint PrepNEMTMC(const uint box, const uint midx = 0, const uint kidx = 0);
  virtual void CalcEn();
  virtual uint Transform();
  virtual void Accept(const uint rejectState, const ulong step);
  virtual void PrintAcceptKind();
*/
private:
  uint bPick;
  bool initMol;
  //SystemPotential sysPotNew;
  //XYZArray molTorqueRef;
  //XYZArray molTorqueNew;
  //XYZArray atomForceRecNew;
  //XYZArray molForceRecNew;
  //XYZArray t_k;
  //XYZArray r_k;
  //Coordinates newMolsPos;
  //COM newCOMs;
  int moveType;
  std::vector<uint> moleculeIndex;
  //const MoleculeLookup& molLookup;
  //Random123Wrapper &r123Wrapper;
  bool allTranslate;

  //VariablesCUDA *cudaVars;
  //CellListGPU *cellListGPU;
  bool isOrthogonal;
  int *kill; // kill the simulation if we started with bad configuration

  double GetCoeff();
  void CalculateTrialDistRot();
  void RotateForceBiased(uint molIndex);
  void TranslateForceBiased(uint molIndex);
  void SetMolInBox(uint box);
  XYZ CalcRandomTransform(XYZ const &lb, double const max, uint molIndex);
  double CalculateWRatio(XYZ const &lb_new, XYZ const &lb_old, XYZ const &k,
                         double max4);
};

inline MultiParticleBrownianGPU::MultiParticleBrownianGPU(System &sys, StaticVals const &statV) :
  MultiParticleBrownian(sys, statV){}

#endif