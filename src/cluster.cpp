/*
 * cluster.cpp
 *
 *  Created on: 12-mar-2018
 *      Author: M. El-Kebir
 */

#include "cluster.h"
#include "dkm/dkm.hpp"

Cluster::Cluster(const Matrix& D,
                 int l)
  : _D(D)
  , _l(l)
  , _mapping(D.getNrCharacters(), -1)
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
  
  auto z = std::get<1>(dkm::kmeans_lloyd(data, _l, seed));
  for (int c = 0; c < n; ++c)
  {
    _mapping[c] = z[c];
  }
}
