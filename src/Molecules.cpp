/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "Molecules.h"
#include "Setup.h"
#include "PDBSetup.h" //For init
#include "MolSetup.h" // ""
#include "FFSetup.h"  // ""
#include "Forcefield.h"
#include "VectorLib.h" //for transfer.
#include <algorithm> //For count.
#include <string>
#include "System.h"

class System;


Molecules::Molecules() : start(NULL), kIndex(NULL), countByKind(NULL),
  chain(NULL), kinds(NULL), pairEnCorrections(NULL),
  pairVirCorrections(NULL), printFlag(true) {}

Molecules::~Molecules(void)
{
  delete[] start;
  delete[] kIndex;
  delete[] countByKind;
  delete[] chain;
  delete[] kinds;
  delete[] pairEnCorrections;
  delete[] pairVirCorrections;
}

void Molecules::Init(Setup & setup, Forcefield & forcefield,
                     System & sys)
{
  pdb_setup::Atoms& atoms = setup.pdb.atoms;
  //Molecule kinds arrays/data.
  count = 0;
  kindsCount = setup.mol.kindMap.size();
  countByKind = new uint[kindsCount];
  kinds = new MoleculeKind[kindsCount];

  // A intermediate vector to store these values so I add kIndex in sorted order.
  std::vector<unsigned int> firstAtomIndex;
  // A intermediate vector to store these values so I add kIndex in sorted order.
  std::vector<unsigned int> kIndexVector;

  uint mk = 0;
  typedef std::map<std::__cxx11::string, mol_setup::MolKind> MolMap;
  for (MolMap::iterator it = setup.mol.kindMap.begin(); it != setup.mol.kindMap.end(); ++it) {

    countByKind[mk] = it->second.kindCount;
    count += it->second.kindCount;

    /* This needs to handle being initialized from frag_x */
    kinds[mk].Init(it->first, setup, forcefield, sys);

    // Insert all of a given molecule kind's instance's firstAtomID's into a sorted vector
    //  It must be sorted so it lines up the xyz data from the pdb file.
    for (std::vector<uint>::iterator nestedIt = it->second.firstAtomID.begin(); 
                      nestedIt != it->second.firstAtomID.end(); ++nestedIt) {
      // Find sorted position
      auto sortedIt = std::upper_bound( firstAtomIndex.begin(), 
                                  firstAtomIndex.end(), *nestedIt); //1
      // Save sorted position
      int sortedPosition = std::distance(sortedIt, firstAtomIndex.begin());
      // Add in sorted position
      firstAtomIndex.insert(sortedIt, *nestedIt); //2 
      // Use the sorted position in firstAtomIndex as a position to add kIndex.
      kIndexVector.insert(kIndexVector.begin() + sortedPosition, mk);
    }
    /* This needs to handle being initialized from frag_x */
    // OLD WAY kinds[mk].Init(atoms.resKindNames[mk]
  }

  if (count == 0) {
    std::cerr << "Error: No Molecules were found in the PSF file(s)!" << std::endl;
    exit(EXIT_FAILURE);
  }

  /* This reflects increasing ordered starting Atom indices ordered to match PDB file(s)*/
  start = new uint [count + 1];
  kIndex = vect::transfer<uint>(kIndexVector);
  start = vect::transfer<uint>(firstAtomIndex);
  chain = vect::transfer<char>(atoms.chainLetter);
  start[count] = atoms.x.size();
  /* This reflect increasing ordered atom indices ordered to match PDB file(s)*/


#if ENSEMBLE == GCMC
  //check to have all the molecules in psf file that defined in config file
  std::map<std::string, double>::const_iterator kindCPIt =
    setup.config.sys.chemPot.cp.begin(),
    lastOne = setup.config.sys.chemPot.cp.end();

  while(kindCPIt != lastOne) {
    std::string molName = kindCPIt->first;
    mol_setup::MolMap::const_iterator dataIterator =
      setup.mol.kindMap.find(molName);
    if(dataIterator == setup.mol.kindMap.end()) {
      std::cerr << "================================================"
                << std::endl << "Error: Molecule " << molName
                << " was not defined in the PDB file." << std::endl
                << "================================================"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    kindCPIt++;
  }
#endif

  if(printFlag) {
    //calculating netcharge of all molecule kind
    double netCharge = 0.0;
    bool hasCharge;
    for (uint mk = 0 ; mk < kindsCount; mk++) {
      netCharge += (countByKind[mk] * kinds[mk].GetMoleculeCharge());
      if(kinds[mk].MoleculeHasCharge()) {
        if(!forcefield.ewald && !forcefield.isMartini) {
          std::cout << "Warning: Charge detected in " << kinds[mk].name
                    << " but Ewald Summaion method is disabled!\n\n";
        } else if(!forcefield.electrostatic && forcefield.isMartini) {
          std::cout << "Warning: Charge detected in " << kinds[mk].name
                    << " but Electrostatic energy calculation is disabled!\n\n";
        }
      }
    }

    if(std::abs(netCharge) > 10E-7) {
      std::cout << "================================================"
                << std::endl << std::endl
                << "Warning: Sum of the charge in the system is: "
                << netCharge << std::endl << std::endl
                <<  "================================================"
                << std::endl;

    }
  }

  //Print LJ nonbonded info
  std::vector<uint> totAtomKind;
  std::vector<std::string> atomNames;
  for(uint i = 0; i < kindsCount; ++i) {
    for(uint j = 0; j < kinds[i].NumAtoms(); j++) {
      if(std::find(totAtomKind.begin(), totAtomKind.end(), kinds[i].AtomKind(j))
          == totAtomKind.end()) {
        totAtomKind.push_back(kinds[i].AtomKind(j));
        atomNames.push_back(kinds[i].atomTypeNames[j]);
      }
    }
  }

  PrintLJInfo(totAtomKind, atomNames, forcefield);
  printFlag = false;

  //Pair Correction matrixes
  pairEnCorrections = new double[kindsCount * kindsCount];
  pairVirCorrections = new double[kindsCount * kindsCount];
  for(uint i = 0; i < kindsCount; ++i) {
    for(uint j = i; j < kindsCount; ++j) {
      pairEnCorrections[i * kindsCount + j] = 0.0;
      pairVirCorrections[i * kindsCount + j] = 0.0;
      for(uint pI = 0; pI < kinds[i].NumAtoms(); ++pI) {
        for(uint pJ = 0; pJ < kinds[j].NumAtoms(); ++pJ) {
          pairEnCorrections[i * kindsCount + j] +=
            forcefield.particles->EnergyLRC(kinds[i].AtomKind(pI),
                                            kinds[j].AtomKind(pJ));
          pairVirCorrections[i * kindsCount + j] +=
            forcefield.particles->VirialLRC(kinds[i].AtomKind(pI),
                                            kinds[j].AtomKind(pJ));
        }
      }
      //set other side of the diagonal
      pairEnCorrections[j * kindsCount + i] =
        pairEnCorrections[i * kindsCount + j];
      pairVirCorrections[j * kindsCount + i] =
        pairVirCorrections[i * kindsCount + j];
    }
  }
}

void Molecules::PrintLJInfo(std::vector<uint> &totAtomKind,
                            std::vector<std::string> &names,
                            Forcefield & forcefield)
{
  if(printFlag) {
    uint size =  totAtomKind.size();
    printf("NonBonded 1-4 parameters:\n");
    if(forcefield.exp6) {
      printf("%-6s %-10s %17s %11s %11s %11s %7s \n", "Type1", "Type2",
             "Epsilon(K)", "Sigma(A)", "Rmax(A)", "Rmin(A)", "alpha");
    } else {
      printf("%-6s %-10s %17s %11s %7s \n", "Type1", "Type2", "Epsilon(K)",
             "Sigma(A)", "N");
    }
    for(uint i = 0; i < size; i++) {
      for(uint j = i; j < size; j++) {
        if(forcefield.exp6) {
          printf("%-6s %-10s %17.4f %11.4f %11.4f %11.4f %7.2f \n",
                 names[i].c_str(), names[j].c_str(),
                 forcefield.particles->GetEpsilon_1_4(totAtomKind[i],
                     totAtomKind[j]),
                 forcefield.particles->GetSigma_1_4(totAtomKind[i],
                     totAtomKind[j]),
                 forcefield.particles->GetRmax_1_4(totAtomKind[i],
                     totAtomKind[j]),
                 forcefield.particles->GetRmin_1_4(totAtomKind[i],
                     totAtomKind[j]),
                 forcefield.particles->GetN_1_4(totAtomKind[i],
                                                totAtomKind[j]));
        } else {
          printf("%-6s %-10s %17.4f %11.4f %7.2f \n", names[i].c_str(),
                 names[j].c_str(),
                 forcefield.particles->GetEpsilon_1_4(totAtomKind[i],
                     totAtomKind[j]),
                 forcefield.particles->GetSigma_1_4(totAtomKind[i],
                     totAtomKind[j]),
                 forcefield.particles->GetN_1_4(totAtomKind[i],
                                                totAtomKind[j]));
        }
      }
    }

    std::cout << std::endl;
    printf("NonBonded parameters:\n");
    if(forcefield.exp6) {
      printf("%-6s %-10s %17s %11s %11s %11s %7s \n", "Type1", "Type2",
             "Epsilon(K)", "Sigma(A)", "Rmax(A)", "Rmin(A)", "alpha");
    } else {
      printf("%-6s %-10s %17s %11s %7s \n", "Type1", "Type2", "Epsilon(K)",
             "Sigma(A)", "N");
    }
    for(uint i = 0; i < size; i++) {
      for(uint j = i; j < size; j++) {
        if(forcefield.exp6) {
          printf("%-6s %-10s %17.4f %11.4f %11.4f %11.4f %7.2f \n",
                 names[i].c_str(), names[j].c_str(),
                 forcefield.particles->GetEpsilon(totAtomKind[i],
                                                  totAtomKind[j]),
                 forcefield.particles->GetSigma(totAtomKind[i],
                                                totAtomKind[j]),
                 forcefield.particles->GetRmax(totAtomKind[i],
                                               totAtomKind[j]),
                 forcefield.particles->GetRmin(totAtomKind[i],
                                               totAtomKind[j]),
                 forcefield.particles->GetN(totAtomKind[i],
                                            totAtomKind[j]));
        } else {
          printf("%-6s %-10s %17.4f %11.4f %7.2f \n", names[i].c_str(),
                 names[j].c_str(),
                 forcefield.particles->GetEpsilon(totAtomKind[i], totAtomKind[j]),
                 forcefield.particles->GetSigma(totAtomKind[i], totAtomKind[j]),
                 forcefield.particles->GetN(totAtomKind[i], totAtomKind[j]));
        }
      }
    }
    std::cout << std::endl;
  }
}
