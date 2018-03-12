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
  /// @param n Number of SNV clusters
  /// @param seed Random number generator seed
  CoordinateAscent(const Matrix& D,
                   int k,
                   bool lazy,
                   double alpha,
                   double beta,
                   int l,
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
  
  const StlIntVector& getZ() const
  {
    return _z;
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
  
  double solveZ();
  
  double computeLogLikelihood(int c, int f) const;
  
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
  /// Number of SNV clusters
  const int _l;
  /// Random number generator seed
  const int _seed;
  /// Output matrix
  Matrix _E;
  /// Character to cluster assignment;
  StlIntVector _z;
  /// Log likelihood
  double _L;
};

#endif // COORDINATEASCENT_H
