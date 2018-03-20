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
  Comparison(const PhylogeneticTree& trueT,
             const PhylogeneticTree& inferredT);
  
  int getRF() const
  {
    return _symDiffSplits.size();
  }
  
  double getLogLikelihood(double alpha, double beta) const;
  
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
  
  void computeLossPrecisionAndRecall(double& precision,
                                     double& recall,
                                     double& lossF1) const;
  
  void computeFlips(const Matrix& input,
                    int& flip01_correct, int& flip01_incorrect,
                    int& flip10_correct, int& flip10_incorrect) const;
  
private:
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
