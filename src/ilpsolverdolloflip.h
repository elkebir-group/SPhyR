/*
 * ilpsolverdolloflip.h
 *
 *  Created on: 25-feb-2018
 *      Author: M. El-Kebir
 */

#ifndef ILPSOLVERDOLLOFLIP_H
#define ILPSOLVERDOLLOFLIP_H

#include "ilpsolverdollo.h"

/// Dollo phylogeny with errors
class IlpSolverDolloFlip : public IlpSolverDollo
{
public:
  /// Constructor
  ///
  /// @param D Input matrix
  /// @param k Maximum number of losses
  /// @param alpha False positve rate (substitution error)
  /// @param beta False negative rate (allele dropout)
  IlpSolverDolloFlip(const Matrix& D,
                     int k,
                     double alpha,
                     double beta);
  
  /// Initialize solver
  virtual void init();
  
  /// Write solution
  ///
  /// @param out Output stream
//  virtual void printSolution(std::ostream& out) const;
  
protected:
  /// Constructor
  ///
  /// @param D Input matrix
  /// @param k Maximum number of losses
  /// @param n Number of characters
  /// @param alpha False negative rate
  /// @param beta False positive rate
  IlpSolverDolloFlip(const Matrix& D,
                     int k,
                     int n,
                     double alpha,
                     double beta);
  
  /// Initialize objective function (number of losses)
  virtual void initObjective();
  
  /// Generate lambda
  ///
  /// @param alpha False positve rate (substitution error)
  /// @param beta False negative rate (allele dropout)
  static double getLambda(double alpha, double beta)
  {
    assert(alpha != 0 || beta != 0);
    assert(0 <= alpha && alpha < 0.5);
    assert(0 <= beta  && beta  < 0.5);
    
    if (alpha == 0)
    {
      // lim_{alpha -> 0} g(alpha,beta) = 1
      return 1;
    }
    else if (beta == 0)
    {
      // lim_{beta -> 0} g(alpha,beta) = 0
      return 0;
    }
    else
    {
      return log( (1 - beta) / alpha ) / log( ((1 - alpha) * (1 - beta)) / (alpha * beta) );
    }
  }
  
protected:
  /// False positve rate (substitution error)
  const double _alpha;
  /// False negative rate (allelic dropout)
  const double _beta;
};

#endif // ILPSOLVERDOLLOFLIP_H
