/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef BONDADJACENCYLIST_H
#define BONDADJACENCYLIST_H

#include <vector>
#include <cstdio>
#include "BasicTypes.h"
#include <iostream>
#include <limits.h>


class Vertex;

class BondAdjacencyList{
public:
BondAdjacencyList(FILE* psf, uint nAtoms, uint nBonds);
//~BondAdjacencyList();

std::vector< std::vector<int> > adjList;
std::vector<Vertex> moleculeConnectedGraph;

private:
int buildAdjList(FILE* psf, uint nBonds);
int buildConnectedGraph(uint nAtoms);
int molCount;
};

class BondAdjacencyMatrix{
public:
BondAdjacencyMatrix(FILE* psf, uint nAtoms, uint nBonds);
//~BondAdjacencyList();

int ** adjMatrix;
int ** degMatrix;
int ** laplacianMatrix;


private:
int buildAdjMatrix(FILE* psf, uint nBonds, uint nAtoms);
int buildDegMatrix(uint nAtoms);
int buildLaplacianMatrix(uint nAtoms);

int molCount;
};

class Vertex{
    public:
Vertex();
//~Vertex();
void Connect(int molCount);
int molID;
bool isConnected;

};
#endif