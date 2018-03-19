/*
 * dollophylogenetictree.h
 *
 *  Created on: 16-mar-2018
 *      Author: M. El-Kebir
 */

#ifndef DOLLOPHYLOGENETICTREE_H
#define DOLLOPHYLOGENETICTREE_H

#include "phylogenetictree.h"
#include "matrix.h"
#include <random>

class DolloPhylogeneticTree : public PhylogeneticTree
{
public:
  DolloPhylogeneticTree(const Matrix& A);
  
  static DolloPhylogeneticTree* parse(const std::string& filename);
  
  bool reconstructTree();
  
  void generateLosses(double lossRate, int k, int seed);
  
  virtual Matrix getMatrixA() const
  {
    return _A;
  }
  
private:
  void expand();
  
  void labelNodes(Node v);
  
  int getOriginalCharacter(int d) const
  {
    const int k = _A.getMaxNrLosses();
    return d / (k + 1);
  }
  
  int getOriginalState(int d) const
  {
    const int k = _A.getMaxNrLosses();
    return (d % (k + 1)) + 1;
  }
  
  int getExpandedCharacter(int c, int i) const
  {
    const int k = _A.getMaxNrLosses();
    assert(1 <= i && i <= k + 1);
    
    return c * (k + 1) + i - 1;
  }
  
  void generateLosses(Node v,
                      int k,
                      double lossRate,
                      StlIntVector& allowedLosses,
                      std::mt19937& rng);
  
private:
  /// Completed matrix
  Matrix _A;
  /// k-Dollo matrix
  Matrix _B;
  /// Binary expansion of _A
  Matrix _Bprime;
  /// Node label (state vectors)
  IntVectorNodeMap _a;
};

#endif // PHYLOGENETICTREE_H

