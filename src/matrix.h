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

/// This class models a (k-Dollo) completion matrix
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
  
  /// Construct matrix from file. Returns NULL if construction fails.
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
  
  /// Return the number of entries with the given value
  ///
  /// @param value Value
  int getCount(int value) const
  {
    int count = 0;
    for (int p = 0; p < _m; ++p)
    {
      for (int c = 0; c < _n; ++c)
      {
        if (getEntry(p, c) == value)
        {
          ++count;
        }
      }
    }
    return count;
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
  
  /// Perturb the entries of the matrix
  ///
  /// @param alpha Probability of changing 0 to 1 (false positive)
  /// @param beta Probability of changing 1 to 0 (false negative)
  /// @param seed Random number generator seed
  void perturb(double alpha, double beta, int seed);
  
  /// Forbidden submatrix
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
    
    /// Character
    const int _c;
    /// Character
    const int _d;
    /// Taxon
    const int _p;
    /// Taxon
    const int _q;
    /// Taxon
    const int _r;
    /// Condition
    const int _condition;
  };
  
  /// List of forbidden submatrices
  typedef std::list<Violation> ViolationList;
  
  /// Identfies forbidden submatrices
  ///
  /// @param int k Maximum number of character losses
  /// @param violationList Output list of forbidden submatrices
  void identifyViolations(int k,
                          ViolationList& violationList) const;
  
  /// Identifies repeated characters (columns)
  ///
  /// @param characterMapping Cluster assignment of original characters
  void identifyRepeatedColumns(StlIntVector& characterMapping) const;
  
  /// Identifies repeated taxa (rows)
  ///
  /// @param taxonMapping Cluster assignment of original taxa
  void identifyRepeatedRows(StlIntVector& taxonMapping) const;
  
  /// Return new matrix with removed repeated and redundant characters and taxa
  ///
  /// @param characterMapping Cluster assignment of original characters
  /// @param taxonMapping Cluster assignment of original taxa
  Matrix simplify(StlIntVector& characterMapping,
                  StlIntVector& taxonMapping) const;
  
  /// Return new matrix with removed repeated and redundant characters
  ///
  /// @param characterMapping Cluster assignment of original characters
  Matrix simplifyColumns(StlIntVector& characterMapping) const;
  
  /// Return new matrix with removed repeated and redundant taxa
  ///
  /// @param taxonMapping Cluster assignment of original taxa
  Matrix simplifyRows(StlIntVector& taxonMapping) const;
  
  /// Return new matrix with previously removed repeated and redundant characters and taxa
  ///
  /// @param characterMapping Cluster assignment of original characters
  /// @param taxonMapping Cluster assignment of original taxa
  Matrix expand(const StlIntVector& characterMapping,
                const StlIntVector& taxonMapping) const;
  
  /// Return new matrix with previously removed repeated and redundant characters
  ///
  /// @param characterMapping Cluster assignment of original characters
  Matrix expandColumns(const StlIntVector& characterMapping) const;

  /// Return new matrix with previously removed repeated and redundant taxa
  ///
  /// @param taxonMapping Cluster assignment of original taxa
  Matrix expandRows(const StlIntVector& taxonMapping) const;
  
  /// Return log Pr(trueMatrix | *this, alpha, beta)
  ///
  /// @param trueMatrix Solution matrix
  /// @param alpha False positive probability
  /// @param beta False negative probability
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

/// Write matrix to output stream
///
/// @param out Output stream
/// @param D Matrix
std::ostream& operator<<(std::ostream& out, const Matrix& D);

/// Read matrix from input stream
///
/// @param in Input stream
/// @param D Matrix
std::istream& operator>>(std::istream& in, Matrix& D);

#endif // MATRIX_H
