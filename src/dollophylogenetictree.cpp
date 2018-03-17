/*
 * dollophylogenetictree.cpp
 *
 *  Created on: 16-mar-2018
 *      Author: M. El-Kebir
 */

#include "dollophylogenetictree.h"

DolloPhylogeneticTree::DolloPhylogeneticTree(const Matrix& A)
  : PhylogeneticTree()
  , _A(A)
  , _B(A.getNrTaxa(), A.getNrCharacters())
  , _Bprime(_A.getNrTaxa(), _A.getNrCharacters() * (_A.getMaxNrLosses() + 1))
  , _a(_T, StlIntVector(_A.getNrCharacters(), 0))
{
  const int m = _A.getNrTaxa();
  const int n = _A.getNrCharacters();
  for (int p = 0; p < m; ++p)
  {
    for (int c = 0; c < n; ++c)
    {
      int a_pc = _A.getEntry(p, c);
      if (a_pc >= 2)
      {
        _B.setEntry(p, c, 0);
      }
    }
  }
  
  expand();
}

DolloPhylogeneticTree* DolloPhylogeneticTree::parse(const std::string& filename)
{
  Matrix* pMatrix = Matrix::parse(filename);
  if (!pMatrix)
  {
    return NULL;
  }
  
  DolloPhylogeneticTree* pTree = new DolloPhylogeneticTree(*pMatrix);
  if (!pTree->reconstructTree())
  {
    std::cerr << "Error: provided matrix '" << filename << "' is not a k-Dollo completion for k = " << pMatrix->getMaxNrLosses() << std::endl;
    delete pMatrix;
    return NULL;
  }
  
  delete pMatrix;
  return pTree;
}

void DolloPhylogeneticTree::expand()
{
  const int m = _A.getNrTaxa();
  const int n = _A.getNrCharacters();
  
  for (int p = 0; p < m; ++p)
  {
    for (int c = 0; c < n; ++c)
    {
      int a_pc = _A.getEntry(p, c);
      if (a_pc == 1 || a_pc >= 2)
      {
        _Bprime.setEntry(p, getExpandedCharacter(c, 1), 1);
      }
      if (a_pc >= 2)
      {
        _Bprime.setEntry(p, getExpandedCharacter(c, a_pc), 1);
      }
    }
  }
}

bool DolloPhylogeneticTree::reconstructTree()
{
  typedef std::vector<IntPair> IntPairVector;
  
  const int m = _Bprime.getNrTaxa();
  const int n = _Bprime.getNrCharacters();
  _taxonToLeaf = NodeVector(m, lemon::INVALID);
  
  StlBoolVector introducedChar(n, false);
  
  // 1. sort columns of Bprime by number of 1s
  IntPairVector columns(n, IntPair(-1, -1));
  for (int c = 0; c < n; ++c)
  {
    int nr1s = 0;
    for (int p = 0; p < m; ++p)
    {
      if (_Bprime.getEntry(p, c) == 1)
      {
        ++nr1s;
      }
    }
    columns[c] = IntPair(nr1s, c);
  }
  // Very important to ensure that in the case of ties losses are introduced after gains!!!
  std::sort(columns.begin(), columns.end(),
            [](const IntPair& a, const IntPair& b)
            {
              return a.first > b.first || (a.first == b.first && a.second < b.second);
            });
  
  for (int p = 0; p < m; ++p)
  {
    Node v = _root;
    for (const IntPair& pair : columns)
    {
      int c = pair.second;

      IntPair charState;
      charState.first = getOriginalCharacter(c);
      charState.second = getOriginalState(c);
      
      if (_Bprime.getEntry(p, c) == 1)
      {
        bool found = false;
        
        for (OutArcIt a(_T, v); a != lemon::INVALID; ++a)
        {
          if (_charStateLabeling[a].count(charState) == 1)
          {
            found = true;
            v = _T.target(a);
          }
        }
        
        if (!found)
        {
          Node w = _T.addNode();
          Arc a = _T.addArc(v, w);
          _charStateLabeling[a].insert(charState);
          _charStateVectorLabeling[a].push_back(charState);
          v = w;
          if (introducedChar[c])
          {
            return false;
          }
          else
          {
            introducedChar[c] = true;
          }
        }
      }
    }
    Node v_p = _T.addNode();
    _T.addArc(v, v_p);
    _taxonToLeaf[p] = v_p;
    _leafToTaxon[v_p] = p;
  }
  
  _a[_root] = StlIntVector(_A.getNrCharacters(), 0);
  _b[_root] = StlIntVector(_A.getNrCharacters(), 0);
  
  labelNodes(_root);
  StlIntVector problem = _b[_taxonToLeaf[0]];
  StlIntVector problem2 = _a[_taxonToLeaf[0]];
  
  collapseBranches();
  
  return true;
}

void DolloPhylogeneticTree::labelNodes(Node v)
{
  for (OutArcIt a(_T, v); a != lemon::INVALID; ++a)
  {
    Node w = _T.target(a);
    _a[w] = _a[v];
    _b[w] = _b[v];
    
    for (const IntPair& ci : _charStateLabeling[a])
    {
      _a[w][ci.first] = ci.second;
      if (ci.second >= 2)
      {
        _b[w][ci.first] = 0;
      }
      else
      {
        assert(ci.second == 1);
        _b[w][ci.first] = 1;
      }
    }
    
    labelNodes(w);
  }
}

void DolloPhylogeneticTree::generateLosses(Node v,
                                           int k,
                                           double lossRate,
                                           StlIntVector&  allowedLosses,
                                           std::mt19937& rng)
{
  std::uniform_real_distribution<> unif(0., 1.);
  const int n = _A.getNrCharacters();
  
  for (OutArcIt a(_T, v); a != lemon::INVALID; ++a)
  {
    Node w = _T.target(a);
    if (OutArcIt(_T, w) == lemon::INVALID)
    {
      // no losses on leaves
      continue;
    }

    // pick a character to lose
    StlIntVector characters;
    for (int c = 0; c < n; ++c)
    {
      if (_a[w][c] == 1
          && _charStateLabeling[a].count(IntPair(c, 1)) == 0
          && allowedLosses[c] > 0)
      {
        characters.push_back(c);
      }
    }
    for (int c : characters)
    {
      double r = unif(rng);
      if (r <= lossRate)
      {
        _charStateLabeling[a].insert(IntPair(c, k - allowedLosses[c] + 2));
        --allowedLosses[c];
      }
    }
  }
  
  labelNodes(v);
  
  for (OutArcIt a(_T, v); a != lemon::INVALID; ++a)
  {
    Node w = _T.target(a);
    
    generateLosses(w, k, lossRate, allowedLosses, rng);
  }
}

void DolloPhylogeneticTree::generateLosses(double lossRate,
                                           int k,
                                           int seed)
{
  StlIntVector allowedLosses(_A.getNrCharacters(), k);
  std::mt19937 rng(seed);
  generateLosses(_root, k, lossRate, allowedLosses, rng);
  
  // Update E and B
  const int m = _Bprime.getNrTaxa();
  const int n = _Bprime.getNrCharacters();
  
  for (int p = 0; p < m; ++p)
  {
    Node v_p = _taxonToLeaf[p];
    for (int c = 0; c < n; ++c)
    {
      _A.setEntry(p, c, _a[v_p][c]);
      _B.setEntry(p, c, _b[v_p][c]);
    }
  }
}
