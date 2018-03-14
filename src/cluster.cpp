/*
 * cluster.cpp
 *
 *  Created on: 12-mar-2018
 *      Author: M. El-Kebir
 */

#include "cluster.h"
#include "dkm/dkm.hpp"

Cluster::Cluster(const Matrix& D,
                 int lT,
                 int lC)
  : _D(D)
  , _lT(lT)
  , _zT(D.getNrTaxa(), -1)
  , _lC(lC)
  , _zC(D.getNrCharacters(), -1)
{
}

void Cluster::cluster(int seed)
{
  // kMeans
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  const int MAX_CELLS = 1000;
  
  assert(m <= MAX_CELLS);
  
  std::vector<std::array<double, MAX_CELLS> > data;
  for (int c = 0; c < n; ++c)
  {
    data.push_back(std::array<double, MAX_CELLS>());
    for (int p = 0; p < m; ++p)
    {
      int d_pc = _D.getEntry(p, c);
      double val = d_pc == -1 ? 0.5 : d_pc;
      data.back()[p] = val;
    }
    for (int p = m; p < MAX_CELLS; ++p)
    {
      data.back()[p] = 0.;
    }
  }
  
  auto z = std::get<1>(dkm::kmeans_lloyd(data, _lC, seed));
  for (int c = 0; c < n; ++c)
  {
    _zC[c] = z[c];
  }
  
  const int MAX_SNVS = 1000;
  
  assert(n <= MAX_SNVS);
  
  std::vector<std::array<double, MAX_SNVS> > data2;
  for (int p = 0; p < m; ++p)
  {
    data2.push_back(std::array<double, MAX_SNVS>());
    for (int c = 0; c < n; ++c)
    {
      int d_pc = _D.getEntry(p, c);
      double val = d_pc == -1 ? 0.5 : d_pc;
      data2.back()[c] = val;
    }
    for (int c = n; c < MAX_SNVS; ++c)
    {
      data2.back()[c] = 0.;
    }
  }
  
  auto z2 = std::get<1>(dkm::kmeans_lloyd(data2, _lT, seed));
  for (int p = 0; p < m; ++p)
  {
    _zT[p] = z2[p];
  }
}
