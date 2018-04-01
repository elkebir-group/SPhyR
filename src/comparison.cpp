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

void Comparison::computeFlips(const Matrix& input,
                              int& flip01_correct, int& flip01_incorrect,
                              int& flip10_correct, int& flip10_incorrect,
                              int& flipx1_correct, int& flipx1_incorrect,
                              int& flipx0_correct, int& flipx0_incorrect) const
{
  Matrix trueB = _trueT.getMatrixB();
  Matrix inferredB = _inferredT.getMatrixB();
  
  assert(trueB.getNrTaxa() == inferredB.getNrTaxa());
  const int m = inferredB.getNrTaxa();
  
  assert(trueB.getNrCharacters() == inferredB.getNrCharacters());
  const int n = inferredB.getNrCharacters();
  
  assert(m == input.getNrTaxa());
  assert(n == input.getNrCharacters());
  
  flip01_correct = flip10_correct = flip01_incorrect = flip10_incorrect = 0;
  flipx1_correct = flipx1_incorrect = flipx0_correct = flipx0_incorrect = 0;
  for (int p = 0; p < m; ++p)
  {
    for (int c = 0; c < n; ++c)
    {
      int input_pc = input.getEntry(p, c);
      int output_pc = inferredB.getEntry(p, c);
      int true_pc = trueB.getEntry(p, c);
      if (input_pc == -1)
      {
        if (output_pc == true_pc)
        {
          if (output_pc == 1)
          {
            ++flipx1_correct;
          }
          else
          {
            ++flipx0_correct;
          }
        }
        else
        {
          if (output_pc == 1)
          {
            ++flipx1_incorrect;
          }
          else
          {
            ++flipx0_incorrect;
          }
        }
      }
      else if (input_pc != output_pc)
      {
        if (output_pc == true_pc)
        {
          if (output_pc == 1)
          {
            ++flip01_correct;
          }
          else
          {
            ++flip10_correct;
          }
        }
        else
        {
          if (output_pc == 1)
          {
            ++flip01_incorrect;
          }
          else
          {
            ++flip10_incorrect;
          }
        }
      }
    }
  }
}

void Comparison::computeLossPrecisionAndRecall(double& precision,
                                               double& recall,
                                               double& lossF1) const
{
  Matrix trueA = _trueT.getMatrixA();
  Matrix inferredA = _inferredT.getMatrixA();
  
  assert(trueA.getNrTaxa() == inferredA.getNrTaxa());
  const int m = inferredA.getNrTaxa();
  
  assert(trueA.getNrCharacters() == inferredA.getNrCharacters());
  const int n = inferredA.getNrCharacters();
  
  int correctlyInferredLosses = 0;
  int inferredLosses = 0;
  int trueLosses = 0;
  
  for (int p = 0; p < m; ++p)
  {
    for (int c = 0; c < n; ++c)
    {
      bool trueLoss = trueA.getEntry(p, c) >= 2;
      if (trueLoss)
      {
        ++trueLosses;
      }
      bool inferredLoss = inferredA.getEntry(p, c) >= 2;
      if (inferredLoss)
      {
        ++inferredLosses;
      }
      if (trueLoss && inferredLoss)
      {
        ++correctlyInferredLosses;
      }
    }
  }
  
  recall = correctlyInferredLosses;
  recall /= trueLosses;
   
  precision = correctlyInferredLosses;
  precision /= inferredLosses;
  
  lossF1 = 2* (precision * recall) / (precision + recall);
}

void Comparison::recallCharStatePairs(double& ancestralRecall,
                                      double& incomparableRecall,
                                      double& clusteredRecall,
                                      bool ignoreLoss) const
{
  PhylogeneticTree::TwoCharStatesSet trueAncestralSet, trueIncomparableSet, trueClusteredSet;
  _trueT.computePairs(trueAncestralSet,
                      trueIncomparableSet,
                      trueClusteredSet,
                      ignoreLoss);
  
  PhylogeneticTree::TwoCharStatesSet inferredAncestralSet, inferredIncomparableSet, inferredClusteredSet;
  _inferredT.computePairs(inferredAncestralSet,
                          inferredIncomparableSet,
                          inferredClusteredSet,
                          ignoreLoss);
  
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
  Matrix trueB = _trueT.getMatrixB();
  Matrix inferredB = _inferredT.getMatrixB();

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
  Matrix trueB = _trueT.getMatrixB();
  Matrix inferredB = _inferredT.getMatrixB();
  
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
