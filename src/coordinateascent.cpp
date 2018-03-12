/*
 * coordinateascent.cpp
 *
 *  Created on: 27-feb-2018
 *      Author: M. El-Kebir
 */

#include "coordinateascent.h"
#include "dkm/dkm.hpp"
#include "ilpsolverdolloflipclustered.h"
#include "columngenflipclustered.h"
#include "cluster.h"

CoordinateAscent::CoordinateAscent(const Matrix& D,
                                   int k,
                                   bool lazy,
                                   double alpha,
                                   double beta,
                                   int l,
                                   int seed)
  : _D(D)
  , _k(k)
  , _lazy(lazy)
  , _alpha(alpha)
  , _beta(beta)
  , _l(l)
  , _seed(seed)
  , _E(_D.getNrTaxa(), l)
  , _z(_D.getNrCharacters(), 0)
  , _L(0)
{
}

void CoordinateAscent::initZ(int seed)
{
  Cluster cluster(_D, _l);
  cluster.cluster(seed);
  _z = cluster.getMapping();
}

double CoordinateAscent::solveE(int timeLimit,
                                int memoryLimit,
                                int nrThreads,
                                bool verbose)
{
//  IlpSolverDolloFlipClustered solvePhylogeny(_D, _k, _alpha, _beta, _l, _z);
  ColumnGenFlipClustered solvePhylogeny(_D, _k, _lazy, _alpha, _beta, _l, _z);
  solvePhylogeny.init();
//  solvePhylogeny.initHotStart(_E, _z);
  solvePhylogeny.solve(timeLimit, memoryLimit, nrThreads, verbose);
  
  _E = solvePhylogeny.getSolA();
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

bool CoordinateAscent::solve(int timeLimit,
                             int memoryLimit,
                             int nrThreads,
                             bool verbose,
                             int nrRestarts)
{
  Matrix bestA;
  double bestLikelihood = -std::numeric_limits<double>::max();
  StlIntVector bestZ;
  
  for (int restartCount = 1; restartCount <= nrRestarts; ++restartCount)
  {
    initZ(_seed + restartCount - 1);
    
    double delta = 1;
    int iteration = 1;
    while (g_tol.nonZero(delta))
    {
      double LL = solveE(timeLimit, memoryLimit, nrThreads, verbose);
      std::cerr << "Restart " << restartCount << " -- iteration " << iteration << " -- E step -- log likelihood " << LL << std::endl;
      double L = solveZ();
      std::cerr << "Restart " << restartCount << " -- iteration " << iteration << " -- z step -- log likelihood " << L << std::endl;
      std::cerr << std::endl;
      
      delta = L - _L;
      _L = L;
      ++iteration;
    }
    
    if (bestLikelihood < _L)
    {
      bestA = _E;
      bestLikelihood = _L;
      bestZ = _z;
    }
  }
  
  _E = bestA;
  _L = bestLikelihood;
  _z = bestZ;
  
  return true;
}
