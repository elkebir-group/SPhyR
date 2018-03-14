/*
 * ilpsolverdolloflipcluster.h
 *
 *  Created on: 27-feb-2018
 *      Author: M. El-Kebir
 */

#ifndef ILPSOLVERDOLLOFLIPCLUSTER_H
#define ILPSOLVERDOLLOFLIPCLUSTER_H

#include "ilpsolverdolloflip.h"

class IlpSolverDolloFlipCluster : public IlpSolverDolloFlip
{
public:
  /// Constructor
  ///
  /// @param D Input matrix
  /// @param k Maximum number of losses
  /// @param alpha False negative rate
  /// @param beta False positive rate
  /// @param lC Number of SNV (characters) clusters
  /// @param lT Number of cell (taxa) clusters
  IlpSolverDolloFlipCluster(const Matrix& D,
                            int k,
                            double alpha,
                            double beta,
                            int lC,
                            int lT);
  
  /// Write solution
  ///
  /// @param out Output stream
  virtual void printSolution(std::ostream& out) const;
  
  /// Initialize cluster assignments
  ///
  /// @param z Cluster assignments
  void initClusters(const StlIntVector& z);
  
  /// Initialize hot start
  ///
  /// @param E Previously identified solution matrix
  /// @param z Previously identified cluster matrix
  void initHotStart(const Matrix& E, const StlIntVector& z);
  
  /// Return cluster assignment
  const StlIntVector& getSolZ() const
  {
    return _solZ;
  }
  
protected:
  /// Initialize variables
  virtual void initVariables();
  
  /// Initialize constraints
  virtual void initConstraints();
  
  /// Initialize objective function
  virtual void initObjective();
  
  /// Extract cplex solution
  virtual void processSolution();
  
protected:
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
  /// Modeling maximum
  IloBoolVarMatrix _w;
  /// Loss
  IloBoolVarMatrix _loss;
  /// Solution cluster assignment
  StlIntVector _solZ;
};

#endif // ILPSOLVERDOLLOFLIPCLUSTER_H
