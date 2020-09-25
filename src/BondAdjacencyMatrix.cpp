/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "BondAdjacencyList.h"


//return if we fail to read anything
const int READERROR = -1;

BondAdjacencyMatrix::BondAdjacencyMatrix(FILE* psf, uint nAtoms, uint nBonds){
    molCount = -1;
    adjMatrix = new int*[nAtoms];
    for (int i = 0; i < nAtoms; i++){
        adjMatrix[i] = new int[nAtoms];
        for (int j = 0; j < nAtoms; j++){
            adjMatrix[i][j] = 0;
        }
    }

    degMatrix = new int*[nAtoms];
    for (int i = 0; i < nAtoms; i++){
        degMatrix[i] = new int[nAtoms];
        for (int j = 0; j < nAtoms; j++){
            degMatrix[i][j] = 0;
        }
    }

    laplacianMatrix = new int*[nAtoms];
    for (int i = 0; i < nAtoms; i++){
        laplacianMatrix[i] = new int[nAtoms];
        for (int j = 0; j < nAtoms; j++){
            laplacianMatrix[i][j] = 0;
        }
    }
    buildAdjMatrix(psf, nBonds, nAtoms);
    buildDegMatrix(nAtoms);
    buildLaplacianMatrix(nAtoms);
}

int BondAdjacencyMatrix::buildAdjMatrix(FILE* psf, uint nBonds, uint nAtoms){
  unsigned int atom0, atom1;
  int dummy;
  for (uint n = 0; n < nBonds; n++) {
    dummy = fscanf(psf, "%u %u", &atom0, &atom1);
    if(dummy != 2) {
      fprintf(stderr, "ERROR: Incorrect Number of bonds in PSF file ");
      return READERROR;
    } else if (feof(psf) || ferror(psf)) {
      fprintf(stderr, "ERROR: Could not find all bonds in PSF file ");
      return READERROR;
    }
    adjMatrix[atom0-1][atom1-1] = 1;
    adjMatrix[atom1-1][atom0-1] = 1;
  }
    std::cout << "\t\tAdjacency Matrix" << std::endl;
    for (int i = 0; i<nAtoms; i++){
        std::cout << "\tcol " << i;
    }
    std::cout << std::endl;
    int row_index = 0;
    for (int i = 0; i < nAtoms; i++){
        std::cout << "row " << i;
        for( int j = 0; j < nAtoms; j++){
            std::cout << "\t" << adjMatrix[i][j];
        }        
        std::cout << std::endl;
    }
}

int BondAdjacencyMatrix::buildDegMatrix(uint nAtoms){
    int degree;
    for (int i = 0; i < nAtoms; i++) {
        degree = 0;
        for (int j = 0; j < nAtoms; j++) {
            degree += adjMatrix[i][j];
        }
        degMatrix[i][i] = degree;
    }
  
    std::cout << "\t\tDegree Matrix" << std::endl;
    for (int i = 0; i<nAtoms; i++){
        std::cout << "\tcol " << i;
    }
    std::cout << std::endl;
    int row_index = 0;
    for (int i = 0; i < nAtoms; i++){
        std::cout << "row " << i;
        for( int j = 0; j < nAtoms; j++){
            std::cout << "\t" << degMatrix[i][j];
        }        
        std::cout << std::endl;
    }
}

int BondAdjacencyMatrix::buildLaplacianMatrix(uint nAtoms){
    for (int i = 0; i < nAtoms; i++) {
        for (int j = 0; j < nAtoms; j++) {
            laplacianMatrix[i][j] = degMatrix[i][j] - adjMatrix[i][j];
        }
    }
  
    std::cout << "\t\tLaplacian Matrix" << std::endl;
    for (int i = 0; i<nAtoms; i++){
        std::cout << "\tcol " << i;
    }
    std::cout << std::endl;
    int row_index = 0;
    for (int i = 0; i < nAtoms; i++){
        std::cout << "row " << i;
        for( int j = 0; j < nAtoms; j++){
            std::cout << "\t" << laplacianMatrix[i][j];
        }        
        std::cout << std::endl;
    }
}
/*
BondAdjacencyList::BondAdjacencyList(FILE* psf, uint nAtoms, uint nBonds){
    adjList.resize(nAtoms + 1);
    moleculeConnectedGraph.resize(nAtoms + 1);
    molCount = -1;
    buildAdjList(psf, nBonds);
    buildConnectedGraph(nAtoms);
}

int BondAdjacencyList::buildAdjList(FILE* psf, uint nBonds){
  unsigned int atom0, atom1;
  int dummy;
  for (uint n = 0; n < nBonds; n++) {
    dummy = fscanf(psf, "%u %u", &atom0, &atom1);
    if(dummy != 2) {
      fprintf(stderr, "ERROR: Incorrect Number of bonds in PSF file ");
      return READERROR;
    } else if (feof(psf) || ferror(psf)) {
      fprintf(stderr, "ERROR: Could not find all bonds in PSF file ");
      return READERROR;
    }
    adjList[atom0].push_back(atom1);
  }
}

int BondAdjacencyList::buildConnectedGraph(uint nAtoms){
  // Evaluate each vertex
  // Start at 1 because psf indexes from 1
    for (uint n = 1; n <= nAtoms; n++) {
    std::cout << "Evaluating atom " << n << std::endl;
    // Evaluate if an incoming edge has connected to me
        if(!moleculeConnectedGraph[n].isConnected){
        // Evaluate if any outgoing edges are connected
            int edgeMolID = -1;
            int min = INT_MAX;
            for(int x = 0; x < adjList[n].size(); x++){
                if(moleculeConnectedGraph[adjList[n][x]].isConnected){
                    min = std::min(min, moleculeConnectedGraph[adjList[n][x]].molID);
                    edgeMolID = min;
                }
            }

            // Neither myself nor any of my outgoing vertices are connected
            // Create a new molecule
            if (edgeMolID == -1){
                molCount++;
                moleculeConnectedGraph[n].Connect(molCount);
                for(int j = 0; j < adjList[n].size(); j++){
                    if(!moleculeConnectedGraph[adjList[n][j]].isConnected){
                        moleculeConnectedGraph[adjList[n][j]].Connect(molCount);
                    }
                }
            // At least one of my outgoing vertices are connected
            // Connect and assign their molID to myself and any other
            // unconnected outgoing vertices.
            } else {
                moleculeConnectedGraph[n].Connect(edgeMolID);
                for(int k = 0; k < adjList[n].size(); k++){
                   // if(!moleculeConnectedGraph[adjList[n][k]].isConnected){
                        moleculeConnectedGraph[adjList[n][k]].Connect(edgeMolID);
                   // }
                }
            }
        } else {
            // I am connected.  I evaluate my outgoing vertices
            // and if any are unconnected, I connect them and 
            // assign my molID to them.
            int min = INT_MAX;
            for(int x = 0; x < adjList[n].size(); x++){
                if(moleculeConnectedGraph[adjList[n][x]].isConnected){
                    min = std::min(min, moleculeConnectedGraph[adjList[n][x]].molID);
                }
            }

            min = std::min(moleculeConnectedGraph[n].molID, min);
            int myID = moleculeConnectedGraph[n].molID;
            for(int i = 1; i <= moleculeConnectedGraph.size(); i++){
                if(moleculeConnectedGraph[i].molID == myID)
                    moleculeConnectedGraph[n].Connect(min);
            }

            for(int i = 0; i < adjList[n].size(); i++){
                if(min != moleculeConnectedGraph[adjList[n][i]].molID){
                    myID = moleculeConnectedGraph[adjList[n][i]].molID;
                    for(int j = 1; j <= moleculeConnectedGraph.size(); i++){
                        if(moleculeConnectedGraph[j].molID == myID)
                            moleculeConnectedGraph[j].Connect(min);
                    }
                }
            }
        }
    }
    for(auto&& x: moleculeConnectedGraph){
        std::cout << x.molID << std::endl;
    }
}

Vertex::Vertex() : molID(-1) , isConnected(false) {}

void Vertex::Connect(int molCount){
    molID = molCount;
    isConnected = true;
}
*/

