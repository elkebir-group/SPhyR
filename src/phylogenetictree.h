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

/// This class models a character-based phylogenetic tree
class PhylogeneticTree
{
public:
  /// Constructor
  PhylogeneticTree();
  
  /// Construct phylogenetic tree from file. Returns NULL if construction fails.
  ///
  /// @param filename Filename
  static PhylogeneticTree* parse(const std::string& filename);
  
  /// Write phylogenetic tree in Graphviz DOT format to output stream
  ///
  /// @param out Output stream
  void writeDOT(std::ostream& out) const;
  
  /// Write phylogenetic tree in Graphviz DOT format to output stream
  ///
  /// @param out Output stream
  /// @param taxonLabel Taxon labels
  /// @param characterLabel Character labels
  void writeDOT(std::ostream& out,
                const StringVector& taxonLabel,
                const StringVector& characterLabel) const;
  
  /// Split, i.e. a bipartition of the leaf set
  typedef std::pair<StlIntSet, StlIntSet> Split;
  
  /// Set of all splits
  typedef std::set<Split> SplitSet;
  
  /// Return split set
  SplitSet getSplitSet() const;
  
  /// Pair of character states
  typedef std::pair<IntPair, IntPair> TwoCharStates;
  
  /// Set of pairs of character states
  typedef std::set<TwoCharStates> TwoCharStatesSet;
  
  /// Classifies pairs of character states as ancestral, incomparable or clustered
  ///
  /// @param ancestral Output set composed of ordered pairs ((c,i),(d,j)) of character states such that (c,i) is ancestral to (d,j)
  /// @param incomparable Output set composed of unordered pairs ((c,i),(d,j)) of character states such that the LCA of (c,i) and (d,j) is the root node
  /// @param clustered Output set composed of unordered pairs ((c,i),(d,j)) of character states such that (c,i) and (d,j) label the same edge
  void computePairs(TwoCharStatesSet& ancestral,
                    TwoCharStatesSet& incomparable,
                    TwoCharStatesSet& clustered,
                    bool ignoreLoss) const;
  
  /// Get number of character states (c,1) that label multiple edges of the tree
  int getParallelEvolutionCount() const;
  
  /// Get number of character states (c,0) that label multiple edges of the tree
  int getBackMutationCount() const;
  
  /// Return binary leaf labeling matrix
  Matrix getMatrixB() const;
  
  /// Return k-Dollo completion of the leaves
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
  
  /// Remove internal nodes with out-degree 1
  void collapseBranches();
  
  /// Compute split set for every node in the subtree rooted at v
  ///
  /// @param v Node
  /// @param universe Character set
  /// @param split Split set
  void computeSplits(Node v,
                     const StlIntSet& universe,
                     SplitNodeMap& split) const;
  
  /// Determine whether a_ci and a_dj are the same
  ///
  /// @param a_ci Edge
  /// @param a_dj Edge
  bool isClustered(Arc a_ci, Arc a_dj) const
  {
    return a_ci == a_dj;
  }
  
  /// Determine whether a_ci is ancestral to a_dj. This is irreflexive.
  ///
  /// @param a_ci Edge
  /// @param a_dj Edge
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
  
  /// Determine whether a_ci and a_dj occur on distinct branches
  ///
  /// @param a_ci Edge
  /// @param a_dj Edge
  bool isIncomparable(Arc a_ci, Arc a_dj) const
  {
    return !isAncestral(a_ci, a_dj) && !isAncestral(a_dj, a_ci);
  }
  
protected:
  /// Phylogenetic tree
  Digraph _T;
  /// Root
  Node _root;
  /// Node label (binary state vectors)
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

/// Write phylogenetic tree to output stream
///
/// @param out Output stream
/// @param T Phylogenetic tree
std::ostream& operator<<(std::ostream& out, const PhylogeneticTree& T);

/// Read phylogenetic tree from input stream
///
/// @param in Input stream
/// @param T Phylogenetic tree
std::istream& operator>>(std::istream& in, PhylogeneticTree& T);

#endif // PHYLOGENETICTREE_H
