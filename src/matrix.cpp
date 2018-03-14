/*
 * matrix.cpp
 *
 *  Created on: 22-feb-2018
 *      Author: M. El-Kebir
 */

#include "matrix.h"
#include <random>

Matrix::Matrix()
  : _m(0)
  , _n(0)
  , _D()
{
}

Matrix::Matrix(int m, int n)
  : _m(m)
  , _n(n)
  , _D(m, StlIntVector(n, 0))
{
}

Matrix Matrix::simplify(StlIntVector& characterMapping,
                        StlIntVector& taxonMapping) const
{
  Matrix newB = simplifyColumns(characterMapping);
  newB = newB.simplifyRows(taxonMapping);
  
  return newB;
}

Matrix Matrix::simplifyRows(StlIntVector& mapping) const
{
  mapping = StlIntVector(_m, -3);
  
  typedef std::set<int> StlIntSet;
  typedef std::pair<int, StlIntSet> IntSetPair;
  typedef std::map<StlIntVector, IntSetPair> StlIntVectorMap;
  
  StlIntVectorMap duplicateMapping;
  
  // identify redundant columns
  for (int p = 0; p < _m; p++)
  {
    int nrZeros = 0;
    int nrOnes = 0;
    int nrMissing = 0;
    
    int lastOne = -1;
    
    for (int c = 0; c < _n; c++)
    {
      int b_pc = getEntry(p, c);
      if (b_pc == 0)
      {
        ++nrZeros;
      }
      else if (b_pc == 1)
      {
        ++nrOnes;
        lastOne = p;
      }
      else
      {
        ++nrMissing;
      }
    }
    
    if (nrZeros + nrMissing == _m)
    {
      mapping[p] = -1; // ALL ZEROS
    }
    else
    {
      StlIntVector vec_p(_n);
      for (int c = 0; c < _n; c++)
      {
        vec_p[c] = getEntry(p, c);
      }
      if (duplicateMapping.count(vec_p) == 0)
      {
        duplicateMapping[vec_p].first = duplicateMapping.size() - 1;
      }
      duplicateMapping[vec_p].second.insert(p);
      mapping[p] = duplicateMapping[vec_p].first;
    }
  }
  
  Matrix newB(duplicateMapping.size(), _n);
  for (int p = 0; p < _m; ++p)
  {
    int pp = mapping[p];
    if (pp >= 0)
    {
      for (int c = 0; c < _n; c++)
      {
        newB.setEntry(pp, c, getEntry(p, c));
      }
    }
  }
  
  return newB;
}

Matrix Matrix::simplifyColumns(StlIntVector& mapping) const
{
  mapping = StlIntVector(_n, -3);
  
  typedef std::set<int> StlIntSet;
  typedef std::pair<int, StlIntSet> IntSetPair;
  typedef std::map<StlIntVector, IntSetPair> StlIntVectorMap;
  
  StlIntVectorMap duplicateMapping;
  
  // identify redundant columns
  for (int c = 0; c < _n; c++)
  {
    int nrZeros = 0;
    int nrOnes = 0;
    int nrMissing = 0;
    
    int lastOne = -1;
    
    for (int p = 0; p < _m; p++)
    {
      int b_pc = getEntry(p, c);
      if (b_pc == 0)
      {
        ++nrZeros;
      }
      else if (b_pc == 1)
      {
        ++nrOnes;
        lastOne = p;
      }
      else
      {
        ++nrMissing;
      }
    }
    
    if (nrZeros + nrMissing == _m)
    {
      mapping[c] = -1; // ALL ZEROS
    }
    else if (nrOnes + nrMissing == _m)
    {
      mapping[c] = -2; // ALL ONES
    }
    else if (nrOnes == 1)
    {
      assert(0 <= lastOne && lastOne < _m);
      mapping[c] = -4 - lastOne; // SINGLE_ONE
    }
    else
    {
      StlIntVector vec_c(_m);
      for (int p = 0; p < _m; p++)
      {
        vec_c[p] = getEntry(p, c);
      }
      if (duplicateMapping.count(vec_c) == 0)
      {
        duplicateMapping[vec_c].first = duplicateMapping.size() - 1;
      }
      duplicateMapping[vec_c].second.insert(c);
      mapping[c] = duplicateMapping[vec_c].first;
    }
  }
  
  Matrix newB(_m, duplicateMapping.size());
  for (int c = 0; c < _n; c++)
  {
    int cc = mapping[c];
    if (cc >= 0)
    {
      for (int p = 0; p < _m; ++p)
      {
        newB.setEntry(p, cc, getEntry(p, c));
      }
    }
  }
  
  return newB;
}


Matrix Matrix::expand(const StlIntVector& characterMapping,
                      const StlIntVector& taxonMapping) const
{
  Matrix newB = expandRows(taxonMapping);
  newB = newB.expandColumns(characterMapping);
  
  return newB;
}

Matrix Matrix::expandRows(const StlIntVector& mapping) const
{
  Matrix newB(mapping.size(), _n);
  
  const int newNrTaxa = mapping.size();
  for (int p = 0; p < newNrTaxa; p++)
  {
    if (mapping[p] == -1)
    {
      for (int c = 0; c < _n; ++c)
      {
        newB.setEntry(p, c, 0);
      }
    }
    else
    {
      for (int c = 0; c < _n; ++c)
      {
        newB.setEntry(p, c, getEntry(mapping[p], c));
      }
    }
  }
  
  return newB;
}

Matrix Matrix::expandColumns(const StlIntVector& mapping) const
{
  Matrix newB(_m, mapping.size());
  
  const int newNrCharacters = mapping.size();
  for (int c = 0; c < newNrCharacters; c++)
  {
    if (mapping[c] == -1)
    {
      for (int p = 0; p < _m; ++p)
      {
        newB.setEntry(p, c, 0);
      }
    }
    else if (mapping[c] == -2)
    {
      for (int p = 0; p < _m; ++p)
      {
        newB.setEntry(p, c, 1);
      }
    }
    else if (mapping[c] <= -4)
    {
      const int taxon = -4 - mapping[c];
      assert(0 <= taxon && taxon < _m);
      for (int p = 0; p < _m; ++p)
      {
        if (p == taxon)
        {
          newB.setEntry(p, c, 1);
        }
        else
        {
          newB.setEntry(p, c, 0);
        }
      }
    }
    else
    {
      for (int p = 0; p < _m; ++p)
      {
        newB.setEntry(p, c, getEntry(p, mapping[c]));
      }
    }
  }
  
  return newB;
}

void Matrix::perturb(double alpha, double beta, int seed)
{
  std::mt19937 rng(seed);
  std::uniform_real_distribution<> unif(0., 1.);
  
  for (int p = 0; p < _m; ++p)
  {
    for (int c = 0; c < _n; ++c)
    {
      if (_D[p][c] == 0)
      {
        if (unif(rng) < alpha)
        {
          _D[p][c] = 1;
        }
      }
      else if (_D[p][c] == 1)
      {
        if (unif(rng) < beta)
        {
          _D[p][c] = 0;
        }
      }
    }
  }
}

void Matrix::identifyViolations(int k,
                                ViolationList& violationList) const
{
  for (int p = 0; p < _m; ++p)
  {
    for (int q = 0; q < _m; ++q)
    {
      if (q == p) continue;
      for (int r = 0; r < _m; ++r)
      {
        if (r == p || r == q) continue;
        for (int c = 0; c < _n; ++c)
        {
          for (int d = c + 1; d < _n; ++d)
          {
            // condition 1
            for (int j = 2; j <= k + 1; ++j)
            {
              for (int j_prime = 1; j_prime <= k + 1; ++j_prime)
              {
                if (j_prime == j) continue;
                
                if (_D[p][c] == 1 && _D[q][c] == 0 && _D[r][c] == 1
                    && _D[p][d] == j_prime && _D[q][d] == j && _D[r][d] == j)
                {
                  violationList.push_back(Violation(c, d, p, q, r, 1));
                }
              }
            }

            // condition 2
            for (int i = 2; i <= k + 1; ++i)
            {
              for (int i_prime = 1; i_prime <= k + 1; ++i_prime)
              {
                if (i == i_prime) continue;
                
                if (_D[p][c] == i && _D[q][c] == i_prime && _D[r][c] == i
                    && _D[p][d] == 0 && _D[q][d] == 1 && _D[r][d] == 1)
                {
                  violationList.push_back(Violation(c, d, p, q, r, 2));
                }
              }
            }
  
            // condition 3
            for (int i = 2; i <= k + 1; ++i)
            {
              for (int i_prime = 1; i_prime <= k + 1; ++i_prime)
              {
                if (i == i_prime) continue;
                for (int j = 2; j <= k + 1; ++j)
                {
                  for (int j_prime = 1; j_prime <= k + 1; ++j_prime)
                  {
                    if (_D[p][c] == i && _D[q][c] == i_prime && _D[r][c] == i
                        && _D[p][d] == j_prime && _D[q][d] == j && _D[r][d] == j)
                    {
                      violationList.push_back(Violation(c, d, p, q, r, 3));
                    }
                  }
                }
              }
            }
            
            // condition 4
            for (int i = 1; i <= k + 1; ++i)
            {
              for (int i_prime = 1; i_prime <= k + 1; ++i_prime)
              {
                for (int j = 1; j <= k + 1; ++j)
                {
                  for (int j_prime = 1; j_prime <= k + 1; ++j_prime)
                  {
                    if (_D[p][c] == i && _D[q][c] == 0 && _D[r][c] == i_prime
                        && _D[p][d] == 0 && _D[q][d] == j && _D[r][d] == j_prime)
                    {
                      violationList.push_back(Violation(c, d, p, q, r, 4));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

std::ostream& operator<<(std::ostream& out, const Matrix& D)
{
  out << D._m << " #taxa" << std::endl;
  out << D._n << " #characters" << std::endl;
  for (int p = 0; p < D._m; ++p)
  {
    bool first = true;
    for (int c = 0; c < D._n; ++c)
    {
      if (first)
        first = false;
      else
        out << " ";
      
      out << D._D[p][c];
    }
    out << std::endl;
  }
  
  return out;
}

std::istream& operator>>(std::istream& in, Matrix& D)
{
  std::string line;
  getline(in, line);
  
  std::stringstream ss(line);

  int m = -1;
  ss >> m;
  if (m < 0)
  {
    throw std::runtime_error(getLineNumber()
                             + "Error: number of taxa should be positive.");
  }
  D._m = m;

  int n = -1;
  getline(in, line);
  ss.clear();
  ss.str(line);
  ss >> n;
  if (n < 0)
  {
    throw std::runtime_error(getLineNumber()
                             + "Error: number of characters should be positive.");
  }
  D._n = n;

  D._D = StlIntMatrix(m, StlIntVector(n, 0));
  for (int p = 0; p < m; ++p)
  {
    StringVector s;
    getline(in, line);
    boost::split(s, line, boost::is_any_of(" \t"));
    if (s.size() < n)
    {
      throw std::runtime_error(getLineNumber()
                               + "Error: insufficient number of characters.");
    }
    
    for (int c = 0; c < n; ++c)
    {
      int i = atoi(s[c].c_str());
      D._D[p][c] = i;
    }
  }
  
  return in;
}
