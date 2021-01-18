
#include <algorithm>
#include <list>
#include <vector>
#include "BondAdjacencyList.h"

template<int N>
class JohnsonsCircuit{

    std::vector<std::list<int>> B;
    BondAdjacencyList AK;
    std::vector<bool> blocked;
    std::vector<int> stack;
    int s;

    bool Circuit(int V);
    void Unblock(int U);    

    JohnsonsCircuit<N>::JohnsonsCircuit(const mol_setup::MolKind & setupKind, std::vector<uint> & bondCount): AK(N, setupKind, bondCount), blocked(N), B(N) {

    }
};