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
  , _s(lT)
  , _zT(D.getNrTaxa(), -1)
  , _t(lC)
  , _zC(D.getNrCharacters(), -1)
{
}

void Cluster::cluster(int seed)
{
  // kMeans
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  const int MAX_CELLS = 1000;
  
  // assert(m <= MAX_CELLS);
  
  // std::vector<std::array<double, MAX_CELLS> > data;
  std::vector<std::vector<double> > data;
  for (int c = 0; c < n; ++c)
  {
    data.push_back(std::vector<double>());
    // data.push_back(std::array<double, MAX_CELLS>());
    for (int i = 0; i < m + 10; ++i)data[c].push_back(0);
    for (int p = 0; p < m; ++p)
    {
      int d_pc = _D.getEntry(p, c);
      double val = d_pc == -1 ? 0.5 : d_pc;
      data.back()[p] = val;
    }
    // for (int p = m; p < MAX_CELLS; ++p)
    // {
    //   data.back()[p] = 0.;
    // }
  }
  
  auto z = std::get<1>(dkm::kmeans_lloyd(data, _t, seed));
  for (int c = 0; c < n; ++c)
  {
    _zC[c] = z[c];
  }
  
  const int MAX_SNVS = 100000;
  
  // assert(n <= MAX_SNVS);
  
  // std::vector<std::array<double, MAX_SNVS> > data2;
  std::vector<std::vector<double> > data2;
  for (int p = 0; p < m; ++p)
  {
    // data2.push_back(std::array<double, MAX_SNVS>());
    data2.push_back(std::vector<double>());
    for (int i = 0; i < n + 10; ++i)data2[p].push_back(0);
    for (int c = 0; c < n; ++c)
    {
      int d_pc = _D.getEntry(p, c);
      double val = d_pc == -1 ? 0.5 : d_pc;
      data2.back()[c] = val;
    }
    // for (int c = n; c < MAX_SNVS; ++c)
    // {
    //   data2.back()[c] = 0.;
    // }
  }
  
  auto z2 = std::get<1>(dkm::kmeans_lloyd(data2, _s, seed));
  for (int p = 0; p < m; ++p)
  {
    _zT[p] = z2[p];
  }
}
