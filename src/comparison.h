/*
 * comparison.h
 *
 *  Created on: 14-mar-2018
 *      Author: M. El-Kebir
 */

#ifndef COMPARISON_H
#define COMPARISON_H

#include "matrix.h"
#include "phylogenetictree.h"

class Comparison
{
public:
  Comparison(const Matrix& trueA,
             const Matrix& inferredA,
             const PhylogeneticTree& trueT,
             const PhylogeneticTree& inferredT);
  
  void compare();
  
  int getTrueK() const
  {
    return _trueA.getMaxNrLosses();
  }
  
  int getInferredK() const
  {
    return _inferredA.getMaxNrLosses();
  }
  
  int getRF() const
  {
    return _symDiffSplits.size();
  }
  
  double getNormalizedRF() const
  {
    PhylogeneticTree::SplitSet combined = _trueSplits;
    combined.insert(_inferredSplits.begin(), _inferredSplits.end());
    return (double) _symDiffSplits.size() / (double) combined.size();
  }
  
  void getTaxaClusteringMetrics(double& RI, double& recall, double& precision) const;
  
  void getCharactersClusteringMetrics(double& RI, double& recall, double& precision) const;
  
  void recallCharStatePairs(double& ancestralRecall,
                            double& incomparableRecall,
                            double& clusteredRecall) const;
  
private:
  /// True matrix
  const Matrix& _trueA;
  /// Inferred matrix
  const Matrix& _inferredA;
  /// True phylogenetic tree
  const PhylogeneticTree& _trueT;
  /// Inferred phylogenetic tree
  const PhylogeneticTree& _inferredT;
  /// True splits
  PhylogeneticTree::SplitSet _trueSplits;
  /// Inferred splits
  PhylogeneticTree::SplitSet _inferredSplits;
  /// Symmetric difference splits
  PhylogeneticTree::SplitSet _symDiffSplits;
};

#endif // COMPARISON_H
