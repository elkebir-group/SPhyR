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

/// This class computes statistics of a solution
class Comparison
{
public:
  /// Constructor
  ///
  /// @param trueT True phylogenetic tree
  /// @param inferredT Inferred phylogenetic tree
  Comparison(const PhylogeneticTree& trueT,
             const PhylogeneticTree& inferredT);
  
  /// Return Robinson-Foulds distance
  int getRF() const
  {
    return _symDiffSplits.size();
  }
  
  /// Return log likelihood of inferred matrix
  ///
  /// @param alpha False positive rate
  /// @param beta False negative rate
  double getLogLikelihood(double alpha, double beta) const;
  
  /// Return normalized Robinson-Foulds distance
  double getNormalizedRF() const
  {
    PhylogeneticTree::SplitSet combined = _trueSplits;
    combined.insert(_inferredSplits.begin(), _inferredSplits.end());
    return (double) _symDiffSplits.size() / (double) combined.size();
  }
  
  /// Return taxa clustering statistics
  ///
  /// @param RI Rand index
  /// @param recall Recall of pairs of taxa
  /// @param precision Precision of pairs of taxa
  void getTaxaClusteringMetrics(double& RI, double& recall, double& precision) const;
  
  /// Return character clustering statistics
  ///
  /// @param RI Rand index
  /// @param recall Recall of pairs of characters
  /// @param precision Precision of pairs of characers
  void getCharactersClusteringMetrics(double& RI, double& recall, double& precision) const;
  
  /// Return tree recall statistics
  ///
  /// @param ancestralRecall Recall of ancestral character states
  /// @param incomparableRecall Recall of incomparable character states
  /// @param clusteredRecall Recall of clustered character states
  /// @param ignoreLoss Ignores character loss in recall measures
  void recallCharStatePairs(double& ancestralRecall,
                            double& incomparableRecall,
                            double& clusteredRecall,
                            bool ignoreLoss) const;
  
  /// Return loss statistics
  ///
  /// @param precision Precision of losses in inferred matrix
  /// @param recall Recall of losses in inferred matrix
  /// @param lossF1 F-1 score of losses in inferred matrix
  void computeLossPrecisionAndRecall(double& precision,
                                     double& recall,
                                     double& lossF1) const;
  
  /// Return flip statistics
  ///
  /// @param input Input matrix
  /// @param flip01_correct Number of correct 0-1 flips
  /// @param flip01_incorrect Number of incorrect 0-1 flips
  /// @param flip10_correct Number of correct 1-0 flips
  /// @param flip10_incorrect Number of incorrect 1-0 flips
  /// @param flipx1_correct Number of correct ?-1 flips
  /// @param flipx1_incorrect Number of incorrect ?-1 flips
  /// @param flip?0_correct Number of correct ?-0 flips
  /// @param flip?0_incorrect Number of incorrect ?-0 flips
  void computeFlips(const Matrix& input,
                    int& flip01_correct, int& flip01_incorrect,
                    int& flip10_correct, int& flip10_incorrect,
                    int& flipx1_correct, int& flipx1_incorrect,
                    int& flipx0_correct, int& flipx0_incorrect) const;
  
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
