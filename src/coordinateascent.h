/*
 * coordinateascent.h
 *
 *  Created on: 27-feb-2018
 *      Author: M. El-Kebir
 */

#ifndef COORDINATEASCENT_H
#define COORDINATEASCENT_H

#include "utils.h"
#include "matrix.h"

/// This class provides a coordinate-ascent based approach to the k-DPFC problem
class CoordinateAscent
{
public:
  /// Constructor
  ///
  /// @param D Input matrix
  /// @param characterMapping Cluster assignment of original characters
  /// @param taxonMapping Cluster assignment of orginal taxa
  /// @param k Maximum number of losses
  /// @param lazy Introduce constraints into the lazy constraint pool
  /// @param alpha False negative rate
  /// @param beta False positive rate
  /// @param s Number of taxon clusters
  /// @param t Number of character clusters
  /// @param seed Random number generator seed
  CoordinateAscent(const Matrix& D,
                   const StlIntVector& characterMapping,
                   const StlIntVector& taxonMapping,
                   int k,
                   bool lazy,
                   double alpha,
                   double beta,
                   int s,
                   int t,
                   int seed);
  
  /// Solve
  ///
  /// @param timeLimit Time limit in seconds
  /// @param memoryLimit Memory limit in megabytes
  /// @param nrThreads Number of threads the solver can use
  /// @param verbose Set to true to enable ILP solver output
  /// @param nrRestart Number of restarts
  bool solve(int timeLimit,
             int memoryLimit,
             int nrThreads,
             bool verbose,
             int nrRestarts);
  
  /// Return solution matrix (k-Dollo completion)
  const Matrix& getE() const
  {
    return _E;
  }
  
  /// Return taxon cluster assignment
  const StlIntVector& getZT() const
  {
    return _zT;
  }
  
  /// Return character cluster assignment
  const StlIntVector& getZC() const
  {
    return _zC;
  }
  
  /// Return log likelihood
  double getLogLikelihood() const
  {
    return _L;
  }
  
private:
  /// Initialize clustering of taxa and characters
  void initZ(int seed);
  
  /// Solve the k-DPFC subproblem given taxon and character clustering. Return log likelihood.
  ///
  /// @param timeLimit Time limit in seconds
  /// @param memoryLimit Memory limit in megabytes
  /// @param nrThreads Number of threads the solver can use
  /// @param verbose Set to true to enable ILP solver output
  /// @param success Indicates whether the optimal solution was found
  double solveE(int timeLimit,
                int memoryLimit,
                int nrThreads,
                bool verbose,
                bool& success);
  
  /// Solve the k-DPFC problem given taxon clustering and k-Dollo completion. Return log likelihood.
  double solveZC();

  /// Solve the k-DPFC problem given character clustering and k-Dollo completion. Return log likelihood.
  double solveZT();
  
  /// Compute log likelihood for output entry (h,f) given that (p,c) is the corresponding input entry.
  ///
  /// @param p Input taxon
  /// @param h Output taxon
  /// @param c Input character
  /// @param f Output character
  double computeLogLikelihood(int p, int h,
                              int c, int f) const;
  
  /// Compute log likelihood for output character f given that c is the corresponding input character
  ///
  /// @param c Input character
  /// @param f Output character
  double computeCharacterLogLikelihood(int c, int f) const;
  
  /// Compute log likelihood for output taxon h given that p is the corresponding input taxon
  ///
  /// @param p Input taxon
  /// @param h Output taxon
  double computeTaxonLogLikelihood(int p, int h) const;
  
  /// Compute log likelihood
  double computeLogLikelihood() const;
  
private:
  /// Input matrix
  const Matrix& _D;
  /// Maximum number of losses
  const int _k;
  /// Use lazy constraints
  const bool _lazy;
  /// False negative rate
  const double _alpha;
  /// False positive rate
  const double _beta;
  /// Number of cell clusters
  const int _s;
  /// Number of character clusters
  const int _t;
  /// Random number generator seed
  const int _seed;
  /// Output matrix
  Matrix _E;
  /// Character to cell assignment
  StlIntVector _zT;
  /// Character to cluster assignment
  StlIntVector _zC;
  /// Log likelihood
  double _L;
  /// Number of correct characters
  double _baseL;
  /// Multiplicative matrix
  StlIntMatrix _multiplicities;
  /// Restart count
  int _restart;
};

#endif // COORDINATEASCENT_H
