/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "MoleculeKind.h"
#include "MolSetup.h"
#include "FFSetup.h"
#include "Forcefield.h"
#include "PDBConst.h" //For resname length.
#include "PRNG.h"
#include "Geometry.h"
#include "Setup.h"
#include "CBMC.h"

#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <algorithm>
#include <utility>
#include <cstdio>
#include <cstdlib>      //for exit



void MoleculeKind::Init
(std::string const& l_name, Setup const& setup,
 Forcefield const& forcefield, System& sys)
{
  mol_setup::MolMap::const_iterator dataIterator =
    setup.mol.kindMap.find(l_name);
  if(dataIterator == setup.mol.kindMap.end()) {
    std::cerr << "================================================"
              << std::endl << "Error: Molecule " << l_name
              << " not found in PDB file. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }
  const mol_setup::MolKind& molData = dataIterator->second;
  name = l_name;

#if ENSEMBLE == GCMC
  std::map<std::string, double>::const_iterator kindCPIt =
    setup.config.sys.chemPot.cp.find(name),
    lastOne = setup.config.sys.chemPot.cp.end();

  //If we don't find a chemical potential for a kind in GCMC mode,
  //then quit.
  if (kindCPIt == lastOne) {
    std::cerr << "================================================"
              << std::endl << "Error: chemical potential is missing for "
              << name << "." << std::endl << std::endl
              << "Here are the listed chemical potentials:"
              << std::endl
              << "----------------------------------------"
              << std::endl;

    //Print out whatever chemical potentials were read.
    for (kindCPIt = setup.config.sys.chemPot.cp.begin();
         kindCPIt != lastOne;
         ++kindCPIt) {
      std::cerr << "Resname: " << kindCPIt->first
                << "      Value: " << kindCPIt->second
                << std::endl;
    }
    exit(EXIT_FAILURE);
  } else {
    chemPot = kindCPIt->second;
  }
#endif

  InitAtoms(molData);

  /* 
   Bonds being numbered from 0 .. n is neccessary for building the graph of the molecule.
   We have bonds and atoms numbered by the unique atom ID from the psf file.  We also already built
   a graph of the simulation, and we found the connected components.  There is nothing stopping 
   us from simply taking the heads of the connected components and deep copying them into a subgraph.
   However, the atom & bond indices would be noncontinous and undefined.  Instead of changing the 
   rest of GOMC to use these new indices, we will renumber them from 0 to n.

   Since we can't guaruntee continuity in the bond numbering, i.e. if we have a disulfide bond 
   between far off amino acids, we can't get clever with
   the smallest atom number by subtracting that number from all atom and bond indices.  
   Therefore, we have to iterate through the bonds in an molecule and renumber them one by one.
  */
  bondList.Init(molData.bonds);
  angles.Init(molData.angles, bondList);
  dihedrals.Init(molData.dihedrals, bondList);

  CF = new CircuitFinder(molData.atoms.size());
  bondCount.resize(NumAtoms(), 0);
  // Count the number of bonds for each atom
  // Add edges to molecule graph
  for (uint i = 0; i < molData.bonds.size(); ++i) {
    ++bondCount[bondList.part1[i]];
    ++bondCount[bondList.part2[i]];
    CF->addEdge(bondList.part1[i], bondList.part2[i]);
    CF->addEdge(bondList.part2[i], bondList.part1[i]);
  }
  //Once-through topology objects

  oneThree = oneFour = false;

  if(forcefield.OneThree) {
    nonBonded_1_3.Init(molData);
    oneThree = true;
  }
  if(forcefield.OneFour) {
    nonBonded_1_4.Init(molData);
    oneFour = true;
  }
  if (forcefield.OneN) {
    nonBonded.Init(molData);
  }

  nonEwaldBonded.Init(molData);
  sortedNB_1_3.Init(nonBonded_1_3, numAtoms);
  sortedNB_1_4.Init(nonBonded_1_4, numAtoms);
  sortedNB.Init(nonBonded, numAtoms);
  sortedEwaldNB.Init(nonEwaldBonded, numAtoms);


#ifdef VARIABLE_PARTICLE_NUMBER
  builder = cbmc::MakeCBMC(sys, forcefield, *this, setup);
  //builder = new cbmc::LinearVlugt(sys, forcefield, *this, setup);
#endif
}

MoleculeKind::MoleculeKind() : angles(3), dihedrals(4),
  atomMass(NULL), builder(NULL), atomKind(NULL),
  atomCharge(NULL) {}


MoleculeKind::~MoleculeKind()
{
  delete[] atomKind;
  delete[] atomMass;
  delete[] atomCharge;
  delete builder;
}

void MoleculeKind::InitAtoms(mol_setup::MolKind const& molData)
{
  numAtoms = molData.atoms.size();
  atomKind = new uint[numAtoms];
  atomMass = new double[numAtoms];
  atomCharge = new double[numAtoms];
  molMass = 0;
  atomNames.clear();
  resNames.clear();

  /* These two entries all PSFOutput to 
    correctly assign residueIDs to a map containing
    multi-residue and standard entries.  */
  isMultiResidue = molData.isMultiResidue;
  intraMoleculeResIDs = molData.intraMoleculeResIDs;

  //convert array of structures to structure of arrays
  for(uint i = 0; i < numAtoms; ++i) {
    const mol_setup::Atom& atom = molData.atoms[i];
    atomNames.push_back(atom.name);
    resNames.push_back(atom.residue);
    atomTypeNames.push_back(atom.type);
    atomMass[i] = atom.mass;
    molMass += atom.mass;
    atomCharge[i] = atom.charge;
    atomKind[i] = atom.kind;
  }
}

double MoleculeKind::GetMoleculeCharge()
{
  double netCharge = 0.0;

  for(uint i = 0; i < numAtoms; ++i) {
    //to calculate net charge
    netCharge += atomCharge[i];
  }

  return netCharge;
}

bool MoleculeKind::MoleculeHasCharge()
{
  for(uint i = 0; i < numAtoms; ++i) {
    if(std::abs(atomCharge[i]) > 0.0)
      return true;
  }
  return false;
}
