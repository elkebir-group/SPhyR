/*
 * phylogenetictree.h
 *
 *  Created on: 28-feb-2018
 *      Author: M. El-Kebir
 */

#ifndef PHYLOGENETICTREE_H
#define PHYLOGENETICTREE_H

#include "utils.h"
#include "matrix.h"
#include <random>

class PhylogeneticTree
{
public:
  PhylogeneticTree(const Matrix& A,
                   int k);
  
  bool reconstructTree();
  
  void writeDOT(std::ostream& out) const;
  
  void generateLosses(double lossRate, int k, int seed);
  
  const Matrix& getA() const
  {
    return _A;
  }
  
  const Matrix& getB() const
  {
    return _B;
  }
  
private:
  void expand();
  
  void label(Node v);
  
  int getOriginalCharacter(int d) const
  {
    return d / (_k + 1);
  }
  
  int getOriginalState(int d) const
  {
    return (d % (_k + 1)) + 1;
  }
  
  int getExpandedCharacter(int c, int i) const
  {
    assert(1 <= i && i <= _k + 1);
    
    return c * (_k + 1) + i - 1;
  }
  
  void generateLosses(Node v, int k, double lossRate, StlIntVector& allowedLosses,
                      std::mt19937& rng);
  
  typedef Digraph::ArcMap<IntPairSet> IntPairSetArcMap;
  
  typedef std::vector<Node> NodeVector;
  
  typedef Digraph::NodeMap<StlIntVector> IntVectorNodeMap;
  
private:
  /// Completed matrix
  Matrix _A;
  /// k-Dollo matrix
  Matrix _B;
  /// Maximum number of losses
  int _k;
  /// Binary expansion of _A
  Matrix _Bprime;
  /// Phylogenetic tree
  Digraph _T;
  /// Root
  Node _root;
  /// Node label (state vectors)
  IntVectorNodeMap _a;
  /// Node label (state vectors)
  IntVectorNodeMap _b;
  ///
  IntPairSetArcMap _charStateLabeling;
  ///
  NodeVector _taxonToNode;
};

#endif // PHYLOGENETICTREE_H
