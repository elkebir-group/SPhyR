/*
 * phylogenetictree.cpp
 *
 *  Created on: 28-feb-2018
 *      Author: M. El-Kebir
 */

#include "phylogenetictree.h"

PhylogeneticTree::PhylogeneticTree(const Matrix& A,
                                   int k)
  : _A(A)
  , _B(A)
  , _k(k)
  , _Bprime(_A.getNrTaxa(), _A.getNrCharacters() * (k + 1))
  , _T()
  , _root(_T.addNode())
  , _a(_T, StlIntVector(_A.getNrCharacters(), 0))
  , _b(_T, StlIntVector(_A.getNrCharacters(), 0))
  , _charStateLabeling(_T)
  , _charStateVectorLabeling(_T)
  , _taxonToLeaf(_A.getNrTaxa(), lemon::INVALID)
  , _leafToTaxon(_T, -1)
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

void PhylogeneticTree::expand()
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
  
//  for (int p = 0; p < m; ++p)
//  {
//    for (int d = 0; d < _Bprime.getNrCharacters(); ++d)
//    {
//      std::cerr << _Bprime.getEntry(p, d) << " ";
//    }
//    std::cerr << std::endl;
//  }
}

bool PhylogeneticTree::reconstructTree()
{
  typedef std::vector<IntPair> IntPairVector;
  
  const int m = _Bprime.getNrTaxa();
  const int n = _Bprime.getNrCharacters();
  
  StlBoolVector introducedChar(n, false);
  
  // 1. sort columns of B by number of 1s
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
  std::sort(columns.begin(), columns.end(), std::greater<IntPair>());
  
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
  label(_root);
  
  collapseBranches();
  
  return true;
}

void PhylogeneticTree::label(Node v)
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
    
    label(w);
  }
}

void PhylogeneticTree::generateLosses(Node v,
                                      int k,
                                      double lossRate,
                                      StlIntVector& allowedLosses,
                                      std::mt19937& rng)
{
  std::uniform_real_distribution<> unif(0., 1.);
  const int n = _A.getNrCharacters();
  
  for (OutArcIt a(_T, v); a != lemon::INVALID; ++a)
  {
    Node w = _T.target(a);

    double r = unif(rng);
    if (r <= lossRate)
    {
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
      if (!characters.empty())
      {
        std::uniform_int_distribution<> unif_int(0, characters.size() - 1);
        int c = characters[unif_int(rng)];
        
        _charStateLabeling[a].insert(IntPair(c, k - allowedLosses[c] + 2));
        --allowedLosses[c];
      }
    }
  }
  
  label(v);
  
  for (OutArcIt a(_T, v); a != lemon::INVALID; ++a)
  {
    Node w = _T.target(a);
    
    generateLosses(w, k, lossRate, allowedLosses, rng);
  }
}

void PhylogeneticTree::generateLosses(double lossRate, int k, int seed)
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

void PhylogeneticTree::writeDOT(std::ostream& out) const
{
  out << "digraph T {" << std::endl;
  
  for (NodeIt v(_T); v != lemon::INVALID; ++v)
  {
    out << "\t" << _T.id(v) << " [label=\"\"]" << std::endl;
  }
  
  for (ArcIt a(_T); a != lemon::INVALID; ++a)
  {
    Node u = _T.source(a);
    Node v = _T.target(a);

    out << "\t" << _T.id(u) << " -> " << _T.id(v) << " [label=\"";
    bool first = true;
    for (const IntPair& ci : _charStateLabeling[a])
    {
      if (first)
        first = false;
      else
        out << "\\n";
      out << ci.first << "," << ci.second;
    }
    out << "\"]" << std::endl;
  }
  
//  for (int p = 0; p < m; ++p)
//  {
//    Node v = _taxonToNode[p];
//    out << "\t" << _T.id(v) << " -> " << "taxon_" << p << std::endl;
//  }
  
  out << "}" << std::endl;
}

void PhylogeneticTree::collapseBranches()
{
  bool done = false;
  while (!done)
  {
    done = true;
    for (NodeIt v(_T); v != lemon::INVALID; ++v)
    {
      int outDegree = lemon::countOutArcs(_T, v);
      if (v != _root && outDegree == 1)
      {
        InArcIt parentArc(_T, v);
        Node parentNode = _T.source(parentArc);
        OutArcIt childArc(_T, v);
        Node childNode = _T.target(childArc);
        
        Arc newArc = _T.addArc(parentNode, childNode);
        _charStateLabeling[newArc].insert(_charStateLabeling[parentArc].begin(),
                                          _charStateLabeling[parentArc].end());
        _charStateLabeling[newArc].insert(_charStateLabeling[childArc].begin(),
                                          _charStateLabeling[childArc].end());
        
        _charStateVectorLabeling[newArc].insert(_charStateVectorLabeling[newArc].end(),
                                                _charStateVectorLabeling[parentArc].begin(),
                                                _charStateVectorLabeling[parentArc].end());
        
        _T.erase(v);
        
        done = false;
        break;
      }
    }
  }
}

void PhylogeneticTree::computePairs(TwoCharStatesSet& ancestral,
                                    TwoCharStatesSet& incomparable,
                                    TwoCharStatesSet& clustered) const
{
  const int n = _A.getNrCharacters();
  
  ArcMatrix charStateToArc(n, ArcVector(_k + 2, lemon::INVALID));
  for (ArcIt a(_T); a != lemon::INVALID; ++a)
  {
    for (const IntPair ci : _charStateLabeling[a])
    {
      charStateToArc[ci.first][ci.second] = a;
    }
  }
  
  for (int c = 0; c < n; ++c)
  {
    for (int i = 0; i <= _k + 1; ++i)
    {
      IntPair ci(c, i);
      // don't distinguish between different losses
      IntPair ci_prime = ci;
      if (i > 2)
        ci_prime.second = 2;

      Arc a_ci = charStateToArc[c][i];
      if (a_ci == lemon::INVALID) continue;
      
      for (int d = 0; d < n; ++d)
      {
        for (int j = 0; j <= _k + 1; ++j)
        {
          IntPair dj(d, j);
          // don't distinguish between different losses
          IntPair dj_prime = dj;
          if (j > 2)
            dj_prime.second = 2;

          Arc a_dj = charStateToArc[d][j];
          if (a_dj == lemon::INVALID) continue;
          
          if (isAncestral(a_ci, a_dj))
          {
            ancestral.insert(TwoCharStates(ci_prime, dj_prime));
          }
          else if (isClustered(a_ci, a_dj))
          {
            if (ci < dj)
            {
              clustered.insert(TwoCharStates(ci_prime, dj_prime));
            }
          }
          else if (isIncomparable(a_ci, a_dj))
          {
            incomparable.insert(TwoCharStates(ci_prime, dj_prime));
          }
        }
      }
    }
  }
}

PhylogeneticTree::SplitSet PhylogeneticTree::getSplitSet() const
{
  StlIntSet universe;
  for (int p = 0; p < _A.getNrTaxa(); ++p)
  {
    universe.insert(p);
  }
  
  SplitSet splitSet;
  
  SplitNodeMap splitMap(_T);
  computeSplits(_root, universe, splitMap);
  
  for (NodeIt v(_T); v != lemon::INVALID; ++v)
  {
    if (v == _root) continue;
    if (OutArcIt(_T, v) == lemon::INVALID) continue;
    
    splitSet.insert(splitMap[v]);
  }
  
  return splitSet;
}

void PhylogeneticTree::computeSplits(Node v,
                                     const StlIntSet& universe,
                                     SplitNodeMap& split) const
{
  if (OutArcIt(_T, v) == lemon::INVALID)
  {
    InArcIt a(_T, v);

    const int p = _leafToTaxon[v];

    split[v].first.clear();
    split[v].first.insert(p);
    split[v].second = universe;
    split[v].second.erase(p);
  }
  else
  {
    split[v].first.clear();
    for (OutArcIt a(_T, v); a != lemon::INVALID; ++a)
    {
      Node w = _T.target(a);
      computeSplits(w, universe, split);
      split[v].first.insert(split[w].first.begin(), split[w].first.end());
    }
    std::set_difference(universe.begin(), universe.end(),
                        split[v].first.begin(), split[v].first.end(),
                        std::inserter(split[v].second, split[v].second.begin()));
  }
}
