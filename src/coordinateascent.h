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

class CoordinateAscent
{
public:
  /// Constructor
  ///
  /// @param D Input matrix
  /// @param k Maximum number of losses
  /// @param alpha False negative rate
  /// @param beta False positive rate
  /// @param gamma Loss rate
  /// @param lT Number of cell clusters
  /// @param lC Number of SNV clusters
  /// @param seed Random number generator seed
  CoordinateAscent(const Matrix& D,
                   const StlIntVector& characterMapping,
                   const StlIntVector& taxonMapping,
                   int k,
                   bool lazy,
                   double alpha,
                   double beta,
                   double gamma,
                   int lT,
                   int lC,
                   int seed);
  
  bool solve(int timeLimit,
             int memoryLimit,
             int nrThreads,
             bool verbose,
             int nrRestarts);
  
  const Matrix& getE() const
  {
    return _E;
  }
  
  const StlIntVector& getZT() const
  {
    return _zT;
  }
  
  const StlIntVector& getZC() const
  {
    return _zC;
  }
  
  double getLogLikelihood() const
  {
    return _L;
  }
  
private:
  void initZ(int seed);
  
  double solveE(int timeLimit,
                int memoryLimit,
                int nrThreads,
                bool verbose);
  
  double solveZC();
  
  double solveZT();
  
  double computeLogLikelihood(int p, int h,
                              int c, int f) const;
  
  double computeCharacterLogLikelihood(int c, int f) const;
  
  double computeTaxonLogLikelihood(int p, int h) const;
  
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
  /// Loss rate
  const double _gamma;
  /// Number of cell clusters
  const int _lT;
  /// Number of SNV clusters
  const int _lC;
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
};

#endif // COORDINATEASCENT_H
