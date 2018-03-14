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
  /// @param lC Number of SNV (characters) clusters
  /// @param zC SNV (characters) cluster assignment
  /// @param lT Number of cell (taxa) clusters
  /// @param zT Cell (taxa) cluster assignment
  IlpSolverDolloFlipClustered(const Matrix& D,
                              int k,
                              double alpha,
                              double beta,
                              int lC,
                              const StlIntVector& zC,
                              int lT,
                              const StlIntVector& zT);
  
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
  /// SNV (character) cluster assignment
  const StlIntVector& _zC;
  /// Cell (taxa) cluster assignment
  const StlIntVector& _zT;
  /// Number of SNV (character) clusters
  const double _lC;
  /// Number of cell (taxa) clusters
  const double _lT;
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
