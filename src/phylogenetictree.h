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
  PhylogeneticTree();
  
  static PhylogeneticTree* parse(const std::string& filename);
  
  void writeDOT(std::ostream& out) const;
  
  typedef std::pair<StlIntSet, StlIntSet> Split;
  
  typedef std::set<Split> SplitSet;
  
  SplitSet getSplitSet() const;
  
  typedef std::pair<IntPair, IntPair> TwoCharStates;
  
  typedef std::set<TwoCharStates> TwoCharStatesSet;
  
  void computePairs(TwoCharStatesSet& ancestral,
                    TwoCharStatesSet& incomparable,
                    TwoCharStatesSet& clustered,
                    bool ignoreLoss) const;
  
  int getParallelEvolutionCount() const;
  
  int getBackMutationCount() const;
  
  Matrix getMatrixB() const;
  
  virtual Matrix getMatrixA() const
  {
    return getMatrixB();
  }
  
protected:
  typedef Digraph::ArcMap<IntPairVector> IntPairVectorArcMap;
  
  typedef Digraph::ArcMap<IntPairSet> IntPairSetArcMap;
  
  typedef std::vector<Node> NodeVector;
  
  typedef Digraph::NodeMap<StlIntVector> IntVectorNodeMap;
  
  typedef Digraph::NodeMap<Split> SplitNodeMap;
  
  void collapseBranches();
  
  void computeSplits(Node v,
                     const StlIntSet& universe,
                     SplitNodeMap& split) const;
  
  bool isClustered(Arc a_ci, Arc a_dj) const
  {
    return a_ci == a_dj;
  }
  
  bool isAncestral(Arc a_ci, Arc a_dj) const
  {
    if (a_ci == a_dj) return false;
    
    while (a_ci != a_dj && a_dj != lemon::INVALID)
    {
      Node v = _T.source(a_dj);
      a_dj = InArcIt(_T, v);
    }
    
    return a_ci == a_dj;
  }
  
  bool isIncomparable(Arc a_ci, Arc a_dj) const
  {
    return !isAncestral(a_ci, a_dj) && !isAncestral(a_dj, a_ci);
  }
  
protected:
  /// Phylogenetic tree
  Digraph _T;
  /// Root
  Node _root;
  /// Node label (state vectors)
  IntVectorNodeMap _b;
  /// Character state arc labeling (unordered)
  IntPairSetArcMap _charStateLabeling;
  /// Character state arc labeling (ordered)
  IntPairVectorArcMap _charStateVectorLabeling;
  /// Taxon index to leaf mapping
  NodeVector _taxonToLeaf;
  /// Leaf to taxon index mapping
  IntNodeMap _leafToTaxon;
  
  friend std::ostream& operator<<(std::ostream& out, const PhylogeneticTree& T);
  friend std::istream& operator>>(std::istream& in, PhylogeneticTree& T);
};

std::ostream& operator<<(std::ostream& out, const PhylogeneticTree& T);
std::istream& operator>>(std::istream& in, PhylogeneticTree& T);

#endif // PHYLOGENETICTREE_H
