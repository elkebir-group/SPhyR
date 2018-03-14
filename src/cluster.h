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

class Cluster
{
public:
  /// Cluster matrix D
  Cluster(const Matrix& D,
          int lT,
          int lC);
  
  void cluster(int seed);
  
  const StlIntVector& getCharacterMapping() const
  {
    return _zC;
  }
  
  const StlIntVector& getTaxonMapping() const
  {
    return _zT;
  }

  
private:
  /// Input matrix
  const Matrix& _D;
  /// Number of clusters
  const int _lT;
  /// Cluster assignment
  StlIntVector _zT;
  /// Number of clusters
  const int _lC;
  /// Cluster assignment
  StlIntVector _zC;
};

#endif // CLUSTER_H
