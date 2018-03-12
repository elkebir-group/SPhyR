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
          int l);
  
  void cluster(int seed);
  
  const StlIntVector& getMapping()
  {
    return _mapping;
  }
  
private:
  /// Input matrix
  const Matrix& _D;
  /// Number of clusters
  const int _l;
  /// Cluster assignment
  StlIntVector _mapping;
};

#endif // CLUSTER_H
