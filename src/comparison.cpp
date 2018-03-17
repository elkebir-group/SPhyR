/*
 * comparison.cpp
 *
 *  Created on: 14-mar-2018
 *      Author: M. El-Kebir
 */

#include "comparison.h"
#include "phylogenetictree.h"

Comparison::Comparison(const PhylogeneticTree& trueT,
                       const PhylogeneticTree& inferredT)
  : _trueT(trueT)
  , _inferredT(inferredT)
  , _trueSplits()
  , _inferredSplits()
  , _symDiffSplits()
{
  _trueSplits = _trueT.getSplitSet();
  _inferredSplits = _inferredT.getSplitSet();
  
  std::set_symmetric_difference(_trueSplits.begin(), _trueSplits.end(),
                                _inferredSplits.begin(), _inferredSplits.end(),
                                std::inserter(_symDiffSplits, _symDiffSplits.begin()));
}

void Comparison::recallCharStatePairs(double& ancestralRecall,
                                      double& incomparableRecall,
                                      double& clusteredRecall) const
{
  PhylogeneticTree::TwoCharStatesSet trueAncestralSet, trueIncomparableSet, trueClusteredSet;
  _trueT.computePairs(trueAncestralSet,
                      trueIncomparableSet,
                      trueClusteredSet);
  
  PhylogeneticTree::TwoCharStatesSet inferredAncestralSet, inferredIncomparableSet, inferredClusteredSet;
  _inferredT.computePairs(inferredAncestralSet,
                          inferredIncomparableSet,
                          inferredClusteredSet);
  
  PhylogeneticTree::TwoCharStatesSet commonAncestralSet, commonIncomparableSet, commonClusteredSet;
  std::set_intersection(trueAncestralSet.begin(), trueAncestralSet.end(),
                        inferredAncestralSet.begin(), inferredAncestralSet.end(),
                        std::inserter(commonAncestralSet, commonAncestralSet.begin()));
  
  std::set_intersection(trueIncomparableSet.begin(), trueIncomparableSet.end(),
                        inferredIncomparableSet.begin(), inferredIncomparableSet.end(),
                        std::inserter(commonIncomparableSet, commonIncomparableSet.begin()));
  
  std::set_intersection(trueClusteredSet.begin(), trueClusteredSet.end(),
                        inferredClusteredSet.begin(), inferredClusteredSet.end(),
                        std::inserter(commonClusteredSet, commonClusteredSet.begin()));
  
  ancestralRecall = (double) commonAncestralSet.size() / (double) trueAncestralSet.size();
  incomparableRecall = (double) commonIncomparableSet.size() / (double) trueIncomparableSet.size();
  clusteredRecall = (double) commonClusteredSet.size() / (double) trueClusteredSet.size();
}

void Comparison::getCharactersClusteringMetrics(double& RI,
                                                double& recall,
                                                double& precision) const
{
  Matrix trueB = _trueT.getMatrix();
  Matrix inferredB = _inferredT.getMatrix();

  assert(trueB.getNrCharacters() == inferredB.getNrCharacters());
  const int n = inferredB.getNrCharacters();
  
  StlIntVector trueCharacterClustering, inferredCharacterClustering;
  trueB.identifyRepeatedColumns(trueCharacterClustering);
  inferredB.identifyRepeatedColumns(inferredCharacterClustering);
  
  IntPairSet TP, TN;
  IntPairSet P, N;
  for (int c = 0; c < n; ++c)
  {
    for (int d = c + 1; d < n; ++d)
    {
      if (trueCharacterClustering[c] == trueCharacterClustering[d])
      {
        P.insert(IntPair(c, d));
        if (inferredCharacterClustering[c] == inferredCharacterClustering[d])
        {
          TP.insert(IntPair(c, d));
        }
      }
      else
      {
        N.insert(IntPair(c, d));
        if (inferredCharacterClustering[c] != inferredCharacterClustering[d])
        {
          TN.insert(IntPair(c, d));
        }
      }
    }
  }
    
  RI = TP.size() + TN.size();
  RI /= (n * (n - 1) / 2);

  recall = TP.size() / (double) P.size();
  precision = TN.size() / (double) N.size();
}

void Comparison::getTaxaClusteringMetrics(double& RI,
                                          double& recall,
                                          double& precision) const
{
  Matrix trueB = _trueT.getMatrix();
  Matrix inferredB = _inferredT.getMatrix();
  
  assert(trueB.getNrCharacters() == inferredB.getNrCharacters());
  const int m = trueB.getNrTaxa();
  
  StlIntVector trueTaxaClustering, inferredTaxaClustering;
  trueB.identifyRepeatedRows(trueTaxaClustering);
  inferredB.identifyRepeatedRows(inferredTaxaClustering);
  
  IntPairSet TP, TN;
  IntPairSet P, N;
  for (int p = 0; p < m; ++p)
  {
    for (int q = p + 1; q < m; ++q)
    {
      if (trueTaxaClustering[p] == trueTaxaClustering[q])
      {
        P.insert(IntPair(p,q));
        if (inferredTaxaClustering[p] == inferredTaxaClustering[q])
        {
          TP.insert(IntPair(p, q));
        }
      }
      else
      {
        N.insert(IntPair(p, q));
        if (inferredTaxaClustering[p] != inferredTaxaClustering[q])
        {
          TN.insert(IntPair(p, q));
        }
      }
    }
  }
  
  RI = TP.size() + TN.size();
  RI /= (m * (m - 1) / 2);
  
  recall = TP.size() / (double) P.size();
  precision = TN.size() / (double) N.size();
}
