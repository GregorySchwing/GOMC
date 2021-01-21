//===----------------------------------------------------------------------===//
//
// Written by Xing Mingjie (mingjie.xing@gmail.com).
//
// An implementation of the Johnson's circuit finding algorithm [1].
//
// [1] Donald B. Johnson, Finding all the elementary circuits of a directed
//     graph, SIAM Journal on Computing, 1975.
//
//===----------------------------------------------------------------------===//

#ifndef CIRCUIT_FINDER_H
#define CIRCUIT_FINDER_H

#include <algorithm>
#include <iostream>
#include <list>
#include <vector>

typedef std::list<int> NodeList;


class CircuitFinder
{
  std::vector< std::vector<int> > uniqueCycles;

  std::vector<NodeList> AK;
  std::vector<int> Stack;
  std::vector<bool> Blocked;
  std::vector<NodeList> B;
  int S;

  int V, E;

  std::vector< std::vector<int> > visited;
  std::vector< std::vector<int> > exactly_1_bonds_apart;
  std::vector< std::vector<int> > exactly_2_bonds_apart;
  std::vector< std::vector<int> > exactly_3_bonds_apart;

  void unblock(int U);
  bool circuit(int V);
  void SaveCycle();
  void run();

  void BFSUtil(int V, int next, int depth);

public:
  void addEdge(int src, int dest);

  CircuitFinder(int N)
    : 
    /* Circuit Finding Vars */
    AK(N), Blocked(N), B(N), 
    /* BFS Vars */
    visited(N), exactly_2_bonds_apart(N), exactly_3_bonds_apart(N) {
        V = N;
        E = 0;
  }
  bool haveCommonElements(std::vector<int> first, std::vector<int> second);
  std::vector<int> returnCombinedSet(std::vector<int> first, std::vector<int> second);
  std::vector< std::vector<int> > GetAllCommonCycles();

  void breadthFirstSearch();
  

};

#endif // CIRCUIT_FINDER_H