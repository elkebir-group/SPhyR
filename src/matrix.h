/*
 * matrix.h
 *
 *  Created on: 22-feb-2018
 *      Author: M. El-Kebir
 */

#ifndef MATRIX_H
#define MATRIX_H

#include "utils.h"
#include <list>

class Matrix
{
public:
  /// Constructor
  ///
  /// @param m Number of taxa
  /// @param n Number of characters
  Matrix(int m, int n);
  
  /// Default constructor
  Matrix();
  
  /// Construct matrix from file
  ///
  /// @param filename Filename
  static Matrix* parse(const std::string& filename);
  
  /// Infer confusion matrix
  ///
  /// @param trueMatrix The true solution
  /// @param TN True negative
  /// @param FN False negative
  /// @param FP False positive
  /// @param TP True positive
  void inferConfusionMatrix(const Matrix& trueMatrix,
                            int& TN, int& FN, int& FP, int& TP) const
  {
    assert(trueMatrix._m == _m);
    assert(trueMatrix._n == _n);
    
    TN = FN = FP = TP = 0;
    for (int p = 0; p < _m; ++p)
    {
      for (int c = 0; c < _n; ++c)
      {
        if (getEntry(p, c) == 0)
        {
          // negative
          if (trueMatrix.getEntry(p, c) == 0)
          {
            // true negative
            ++TN;
          }
          else
          {
            // false negative
            ++FN;
          }
        }
        else
        {
          // positive
          if (trueMatrix.getEntry(p, c) == 0)
          {
            // false positive
            ++FP;
          }
          else
          {
            // true positive
            ++TP;
          }
        }
      }
    }
  }
  
  /// Return number of taxa
  int getNrTaxa() const
  {
    return _m;
  }
  
  /// Return number of characters
  int getNrCharacters() const
  {
    return _n;
  }
  
  /// Return maximum number of losses per character
  int getMaxNrLosses() const
  {
    return _k;
  }
  
  /// Return number of ones
  int getNrOfOnes(int c) const
  {
    assert(0 <= c && c < _n);
    
    int res = 0;
    for (int p = 0; p < _m; ++p)
    {
      if (_D[p][c] == 1)
      {
        ++res;
      }
    }
    
    return res;
  }
  
  /// Return input entry
  ///
  /// @param p Taxon
  /// @param c Character
  int getEntry(int p, int c) const
  {
    assert(0 <= p && p < _m);
    assert(0 <= c && c < _n);
    
    return _D[p][c];
  }
  
  /// Set entry
  ///
  /// @param p Taxon
  /// @param c Character
  /// @param i State
  void setEntry(int p, int c, int i)
  {
    assert(0 <= p && p < _m);
    assert(0 <= c && c < _n);

    _D[p][c] = i;
    
    if (i - 1 > _k)
    {
      _k = i - 1;
    }
  }
  
  void perturb(double alpha, double beta, int seed);
  
  struct Violation
  {
  public:
    Violation(int c, int d,
              int p, int q, int r,
              int condition)
      : _c(c)
      , _d(d)
      , _p(p)
      , _q(q)
      , _r(r)
      , _condition(condition)
    {
    }
    
    const int _c;
    const int _d;
    const int _p;
    const int _q;
    const int _r;
    const int _condition;
  };
  
  void identifyRepeatedColumns(StlIntVector& characterMapping) const;
  
  void identifyRepeatedRows(StlIntVector& taxonMapping) const;
  
  typedef std::list<Violation> ViolationList;
  
  void identifyViolations(int k,
                          ViolationList& violationList) const;
  
  Matrix simplify(StlIntVector& characterMapping,
                  StlIntVector& taxonMapping) const;
  
  Matrix simplifyColumns(StlIntVector& characterMapping) const;
  
  Matrix simplifyRows(StlIntVector& taxonMapping) const;
  
  Matrix expand(const StlIntVector& characterMapping,
                const StlIntVector& taxonMapping) const;
  
  Matrix expandColumns(const StlIntVector& characterMapping) const;
  
  Matrix expandRows(const StlIntVector& taxonMapping) const;
  
  double getLogLikelihood(const Matrix& trueMatrix,
                          double alpha,
                          double beta) const;
  
protected:
  /// Number of taxa
  int _m;
  /// Number of characters
  int _n;
  /// Input matrix
  StlIntMatrix _D;
  /// Number of losses
  int _k;

  friend std::ostream& operator<<(std::ostream& out, const Matrix& D);
  friend std::istream& operator>>(std::istream& in, Matrix& D);
};

std::ostream& operator<<(std::ostream& out, const Matrix& D);
std::istream& operator>>(std::istream& in, Matrix& D);

#endif // MATRIX_H
