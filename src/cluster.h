/*
 * cluster.h
 *
 *  Created on: 12-mar-2018
 *      Author: M. El-Kebir
 */

#ifndef CLUSTER_H
#define CLUSTER_H

#include "utils.h"
#include "matrix.h"

/// This class models k-Means clustering of the character and taxon set
class Cluster
{
public:
  /// Cluster characters and taxa of the given matrix
  ///
  /// @param D Input matrix
  /// @param s Number of character clusters
  /// @param t Number of taxon clusters
  Cluster(const Matrix& D,
          int s,
          int t);
  
  /// Cluster
  ///
  /// @param seed Random number generator seed
  void cluster(int seed);
  
  /// Return character cluster assignment
  const StlIntVector& getCharacterMapping() const
  {
    return _zC;
  }
  
  /// Return taxon cluster assignment
  const StlIntVector& getTaxonMapping() const
  {
    return _zT;
  }

private:
  /// Input matrix
  const Matrix& _D;
  /// Number of taxon clusters
  const int _s;
  /// Taxon cluster assignment
  StlIntVector _zT;
  /// Number of character clusters
  const int _t;
  /// Character cluster assignment
  StlIntVector _zC;
};

#endif // CLUSTER_H
