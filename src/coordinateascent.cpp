/*
 * coordinateascent.cpp
 *
 *  Created on: 27-feb-2018
 *      Author: M. El-Kebir
 */

#include "coordinateascent.h"
#include "dkm/dkm.hpp"
//#include "ilpsolverdolloflipclustered.h"
#include "columngenflipclustered.h"
#include "cluster.h"

CoordinateAscent::CoordinateAscent(const Matrix& D,
                                   int k,
                                   bool lazy,
                                   double alpha,
                                   double beta,
                                   int lT,
                                   int lC,
                                   int seed)
  : _D(D)
  , _k(k)
  , _lazy(lazy)
  , _alpha(alpha)
  , _beta(beta)
  , _lT(lT)
  , _lC(lC)
  , _seed(seed)
  , _E(lT, lC)
  , _zT(_D.getNrTaxa(), 0)
  , _zC(_D.getNrCharacters(), 0)
  , _L(0)
{
}

void CoordinateAscent::initZ(int seed)
{
  Cluster cluster(_D, _lT, _lC);
  cluster.cluster(seed);
  _zT = cluster.getTaxonMapping();
  _zC = cluster.getCharacterMapping();
}

double CoordinateAscent::solveE(int timeLimit,
                                int memoryLimit,
                                int nrThreads,
                                bool verbose)
{
//  IlpSolverDolloFlipClustered solvePhylogeny(_D, _k, _alpha, _beta, _l, _z);
  ColumnGenFlipClustered solvePhylogeny(_D, _k, _lazy, _alpha, _beta,
                                        _lC, _zC, _lT, _zT);
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

double CoordinateAscent::computeCharacterLogLikelihood(int c, int f) const
{
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  assert(0 <= c && c < n);
  assert(0 <= f && f < _lC);
  
  double L = 0;
  for (int p = 0; p < m; ++p)
  {
    int h = _zT[p];
    L += computeLogLikelihood(p, h, c, f);
  }
  return L;
}

double CoordinateAscent::computeTaxonLogLikelihood(int p, int h) const
{
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  assert(0 <= p && p < m);
  assert(0 <= h && h < _lT);
  
  double L = 0;
  for (int c = 0; c < n; ++c)
  {
    int f = _zC[c];
    L += computeLogLikelihood(p, h, c, f);
  }
  return L;
}

double CoordinateAscent::computeLogLikelihood(int p, int h,
                                              int c, int f) const
{
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  assert(0 <= c && c < n);
  assert(0 <= f && f < _lC);
  
  assert(0 <= p && p < m);
  assert(0 <= h && h < _lT);
  
  const double log_alpha = log(_alpha);
  const double log_1_minus_alpha = log(1 - _alpha);
  const double log_beta = log(_beta);
  const double log_1_minus_beta = log(1 - _beta);
  
  double L = 0;
  int d_pc = _D.getEntry(p, c);
  int e_hf = _E.getEntry(h, f);
  assert(d_pc == 0 || d_pc == 1 || d_pc == -1);
  
  if (d_pc == 0)
  {
    if (e_hf == 0)
    {
      L += log_1_minus_alpha;
    }
    else if (e_hf == 1)
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
    if (e_hf == 0)
    {
      L += log_alpha;
    }
    else if (e_hf == 1)
    {
      L += log_1_minus_beta;
    }
    else
    {
      L += log_alpha;
    }
  }
  return L;
}

double CoordinateAscent::computeLogLikelihood() const
{
  const int m = _D.getNrTaxa();
  const int n = _D.getNrCharacters();
  
  double L = 0;
  
  for (int p = 0; p < m; ++p)
  {
    const int h = _zT[p];
    for (int c = 0; c < n; ++c)
    {
      const int f = _zC[c];
      L += computeLogLikelihood(p, h, c, f);
    }
  }
  
  return L;
}

double CoordinateAscent::solveZC()
{
  const int n = _D.getNrCharacters();
  
  double L = 0;
  for (int c = 0; c < n; ++c)
  {
    double maxL_c = -std::numeric_limits<double>::max();
    int max_f = -1;
    for (int f = 0; f < _lC; ++f)
    {
      double L_cf = computeCharacterLogLikelihood(c, f);
      if (L_cf > maxL_c)
      {
        maxL_c = L_cf;
        max_f = f;
      }
    }
    
    assert(max_f != -1);
    _zC[c] = max_f;
    L += maxL_c;
  }
  
  return L;
}

double CoordinateAscent::solveZT()
{
  const int m = _D.getNrTaxa();
  
  double L = 0;
  for (int p = 0; p < m; ++p)
  {
    double maxL_p = -std::numeric_limits<double>::max();
    int max_h = -1;
    for (int h = 0; h < _lT; ++h)
    {
      double L_ph = computeTaxonLogLikelihood(p, h);
      if (L_ph > maxL_p)
      {
        maxL_p = L_ph;
        max_h = h;
      }
    }
    
    assert(max_h != -1);
    _zT[p] = max_h;
    L += maxL_p;
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
  StlIntVector bestZT, bestZC;
  
  for (int restartCount = 1; restartCount <= nrRestarts; ++restartCount)
  {
    initZ(_seed + restartCount - 1);
    
    double delta = 1;
    int iteration = 1;
    while (g_tol.nonZero(delta))
    {
      double LLL = solveE(timeLimit, memoryLimit, nrThreads, verbose);
      std::cerr << "Restart " << restartCount << " -- iteration " << iteration << " -- E step -- log likelihood " << LLL << std::endl;
//      std::cout << _E << std::endl;
      
      double LL = solveZT();
      std::cerr << "Restart " << restartCount << " -- iteration " << iteration << " -- zT step -- log likelihood " << LL << std::endl;
      double L = solveZC();
      std::cerr << "Restart " << restartCount << " -- iteration " << iteration << " -- zC step -- log likelihood " << L << std::endl;
      std::cerr << std::endl;
      
      delta = L - _L;
      _L = L;
      ++iteration;
    }
    
    if (bestLikelihood < _L)
    {
      bestA = _E;
      bestLikelihood = _L;
      bestZT = _zT;
      bestZC = _zC;
    }
  }
  
  _E = bestA;
  _L = bestLikelihood;
  _zT = bestZT;
  _zC = bestZC;
  
  return true;
}
