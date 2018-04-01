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

/// This class models a k-Dollo phylogenetic tree
class DolloPhylogeneticTree : public PhylogeneticTree
{
public:
  /// Constructor
  ///
  /// @param A k-Dollo completion
  DolloPhylogeneticTree(const Matrix& A);
  
  /// Construct k-Dollo phylogenetic tree from file. Returns NULL if construction fails.
  ///
  /// @param filename Filename
  static DolloPhylogeneticTree* parse(const std::string& filename);
  
  /// Attempts to construct k-Dollo phylogenetic tree from the input k-Dollo completion
  bool reconstructTree();
  
  /// Introduces character losses
  ///
  /// @param lossRate Probability of character loss
  /// @param k Maximum number of losses per character
  /// @param seed Random number generator seed
  void generateLosses(double lossRate, int k, int seed);
  
  /// Return k-Dollo completion of the leaves
  virtual Matrix getMatrixA() const
  {
    return _A;
  }
  
private:
  /// Generate binary factor matrix
  void expand();
  
  /// Label nodes in the subtree rooted at the given node the character states they posses
  ///
  /// @param v Node
  void labelNodes(Node v);
  
  /// Maps character of binary factor matrix back to original character
  int getOriginalCharacter(int d) const
  {
    const int k = _A.getMaxNrLosses();
    return d / (k + 1);
  }
  
  /// Maps character of binary factor matrix back to original state
  int getOriginalState(int d) const
  {
    const int k = _A.getMaxNrLosses();
    return (d % (k + 1)) + 1;
  }
  
  /// Maps character state of k-Dollo completion to character of binary factor matrix
  ///
  /// @param c Character
  /// @param i State
  int getExpandedCharacter(int c, int i) const
  {
    const int k = _A.getMaxNrLosses();
    assert(1 <= i && i <= k + 1);
    
    return c * (k + 1) + i - 1;
  }
  
  /// Introduce losses at the subtree rooted at the given node
  ///
  /// @param v Node
  /// @param k Maximum number of losses per character
  /// @param lossRate Probability of character loss
  /// @param allowedLosses Number of remaining losses for each character
  /// @param rng Random number generator
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
