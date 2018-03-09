/*
 * ilpsolverdolloflipclustered.h
 *
 *  Created on: 27-feb-2018
 *      Author: M. El-Kebir
 */

#ifndef ILPSOLVERDOLLOFLIPCLUSTERED_H
#define ILPSOLVERDOLLOFLIPCLUSTERED_H

#include "ilpsolverdolloflip.h"

class IlpSolverDolloFlipClustered : public IlpSolverDolloFlip
{
public:
  /// Constructor
  ///
  /// @param D Input matrix
  /// @param k Maximum number of losses
  /// @param alpha False negative rate
  /// @param beta False positive rate
  /// @param l Number of SNV clusters
  /// @param z Cluster assignment
  IlpSolverDolloFlipClustered(const Matrix& D,
                              int k,
                              double alpha,
                              double beta,
                              int l,
                              const StlIntVector& z);
  
  /// Initialize hot start
  ///
  /// @param E Previously identified solution matrix
  /// @param z Previously identified cluster matrix
  void initHotStart(const Matrix& E, const StlIntVector& z);
  
  void initMaxLoss(int maxL);
  
protected:
  /// Initialize variables
  virtual void initVariables();
  
  /// Initialize constraints
  virtual void initConstraints();
  
  /// Initialize objective function
  virtual void initObjective();
  
protected:
  /// Cluster assignment
  const StlIntVector& _z;
  /// Number of SNV clusters
  const double _l;
  /// Best-case likelihoods
  StlDoubleVector _bestCaseLikelihoods;
  /// Worst-case likelihoods
  StlDoubleVector _worstCaseLikelihoods;
  /// Likelihood
  IloNumVarArray _L;
  /// Loss
  IloBoolVarArray _loss;
};

#endif // ILPSOLVERDOLLOFLIPCLUSTERED_H
