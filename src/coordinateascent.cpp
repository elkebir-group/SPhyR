/*
 * coordinateascent.cpp
 *
 *  Created on: 27-feb-2018
 *      Author: M. El-Kebir
 */

#include "coordinateascent.h"
#include "dkm/dkm.hpp"
#include "ilpsolverdolloflipclustered.h"

CoordinateAscent::CoordinateAscent(const Matrix& D,
                                   int k,
                                   double alpha,
                                   double beta,
                                   int l,
                                   int seed)
  : _D(D)
  , _k(k)
  , _alpha(alpha)
  , _beta(beta)
  , _l(l)
  , _seed(seed)
  , _E(_D.getNrTaxa(), l)
  , _z(_D.getNrCharacters(), 0)
  , _L(0)
{
}

void CoordinateAscent::initZ()
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
  
  auto z = std::get<1>(dkm::kmeans_lloyd(data, _l, _seed));
  for (int c = 0; c < n; ++c)
  {
    _z[c] = z[c];
  }
}

double CoordinateAscent::solveE(int timeLimit,
                                int memoryLimit,
                                int nrThreads)
{
  IlpSolverDolloFlipClustered solvePhylogeny(_D, _k, _alpha, _beta, _l, _z);
  solvePhylogeny.init();
  solvePhylogeny.initHotStart(_E, _z);
  solvePhylogeny.solve(timeLimit, memoryLimit, nrThreads);
  
  _E = solvePhylogeny.getSolE();
#ifdef DEBUG
  InputMatrix::ViolationList violationList;
  assert(_E.identifyViolations(_k, violationList));
  assert(violationList.empty());
#endif // DEBUG
  
  return computeLogLikelihood();
}

double CoordinateAscent::computeLogLikelihood(int c, int f) const
{
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  assert(0 <= c && c < n);
  assert(0 <= f && f < _l);
  
  const double log_alpha = log(_alpha);
  const double log_1_minus_alpha = log(1 - _alpha);
  const double log_beta = log(_beta);
  const double log_1_minus_beta = log(1 - _beta);
  
  double L = 0;
  for (int p = 0; p < m; p++)
  {
    int d_pc = _D.getEntry(p, c);
    int e_pf = _E.getEntry(p, f);
    assert(d_pc == 0 || d_pc == 1 || d_pc == -1);
    
    if (d_pc == 0)
    {
      if (e_pf == 0)
      {
        L += log_1_minus_alpha;
      }
      else if (e_pf == 1)
      {
        L += log_beta;
      }
      else
      {
        L += log_1_minus_alpha;
      }
    }
    else if (d_pc == 1)
    {
      if (e_pf == 0)
      {
        L += log_alpha;
      }
      else if (e_pf == 1)
      {
        L += log_1_minus_beta;
      }
      else
      {
        L += log_alpha;
      }
    }
  }
  return L;
}

double CoordinateAscent::computeLogLikelihood() const
{
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  double L = 0;
  
  for (int c = 0; c < n; ++c)
  {
    const int f = _z[c];
    L += computeLogLikelihood(c, f);
  }
  
  return L;
}

double CoordinateAscent::solveZ()
{
  const int n = _D.getNrCharacters();
  
  double L = 0;
  for (int c = 0; c < n; ++c)
  {
    double maxL_c = -std::numeric_limits<double>::max();
    int max_f = -1;
    for (int f = 0; f < _l; ++f)
    {
      double L_cf = computeLogLikelihood(c, f);
      if (L_cf > maxL_c)
      {
        maxL_c = L_cf;
        max_f = f;
      }
    }
    
    assert(max_f != -1);
    _z[c] = max_f;
    L += maxL_c;
  }
  
  return L;
}

double CoordinateAscent::solve(int timeLimit,
                               int memoryLimit,
                               int nrThreads)
{
  initZ();
  
  double delta = 1;
  int iteration = 1;
  while (g_tol.nonZero(delta))
  {
    double LL = solveE(timeLimit, memoryLimit, nrThreads);
    std::cerr << "Iteration " << iteration << " -- E step -- log likelihood " << LL << std::endl;
    double L = solveZ();
    std::cerr << "Iteration " << iteration << " -- z step -- log likelihood " << L << std::endl;
    
    delta = L - _L;
    _L = L;
  }
  
  return _L;
}
