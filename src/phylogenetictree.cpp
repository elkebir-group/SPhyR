/*
 * phylogenetictree.cpp
 *
 *  Created on: 28-feb-2018
 *      Author: M. El-Kebir
 */

#include "phylogenetictree.h"

PhylogeneticTree::PhylogeneticTree()
  : _T()
  , _root(_T.addNode())
  , _b(_T)
  , _charStateLabeling(_T)
  , _charStateVectorLabeling(_T)
  , _taxonToLeaf()
  , _leafToTaxon(_T, -1)
{
}

PhylogeneticTree* PhylogeneticTree::parse(const std::string& filename)
{
  PhylogeneticTree* pTree = NULL;
  
  std::ifstream inT(filename.c_str());
  if (!inT.good())
  {
    std::cerr << "Error: could not open '" << filename << "' for reading" << std::endl;
    return NULL;
  }
  
  pTree = new PhylogeneticTree();
  inT >> *pTree;
  inT.close();
  
  return pTree;
}

void PhylogeneticTree::writeDOT(std::ostream& out) const
{
  out << "digraph T {" << std::endl;
  
  for (NodeIt v(_T); v != lemon::INVALID; ++v)
  {
    if (OutArcIt(_T, v) == lemon::INVALID)
    {
      out << "\t" << _T.id(v) << " [label=\"" << _leafToTaxon[v] << "\"]" << std::endl;
    }
    else
    {
      out << "\t" << _T.id(v) << " [label=\"\"]" << std::endl;
    }
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
  const int n = _b[_root].size();
  
  typedef std::vector<Arc> ArcVector;
  typedef std::vector<ArcVector> ArcMatrix;
  typedef std::vector<ArcMatrix> Arc3Matrix;
  
  const int nrStates = 2;
  
  Arc3Matrix charStateToArc(n, ArcMatrix(nrStates, ArcVector()));
  for (ArcIt a(_T); a != lemon::INVALID; ++a)
  {
    for (IntPair ci : _charStateLabeling[a])
    {
      if (ci.second >= 2)
        ci.second = 0;
      charStateToArc[ci.first][ci.second].push_back(a);
    }
  }
  
  for (int c = 0; c < n; ++c)
  {
    for (int i = 0; i < nrStates; ++i)
    {
      IntPair ci(c, i);

      for (Arc a_ci : charStateToArc[c][i])
      {
        for (int d = 0; d < n; ++d)
        {
          for (int j = 0; j < nrStates; ++j)
          {
            IntPair dj(d, j);
            if (ci == dj) continue;

            for (Arc a_dj : charStateToArc[d][j])
            {
              if (isAncestral(a_ci, a_dj))
              {
                ancestral.insert(TwoCharStates(ci, dj));
              }
              else if (isClustered(a_ci, a_dj))
              {
                if (ci < dj)
                {
                  clustered.insert(TwoCharStates(ci, dj));
                }
              }
              else if (isIncomparable(a_ci, a_dj))
              {
                incomparable.insert(TwoCharStates(ci, dj));
              }
            }
          }
        }
      }
    }
  }
}

PhylogeneticTree::SplitSet PhylogeneticTree::getSplitSet() const
{
  const int m = _taxonToLeaf.size();
  StlIntSet universe;
  for (int p = 0; p < m; ++p)
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

Matrix PhylogeneticTree::getMatrix() const
{
  const int m = _taxonToLeaf.size();
  const int n = _b[_root].size();
  
  Matrix B(m, n);
  for (int p = 0; p < m; ++p)
  {
    Node v_p = _taxonToLeaf[p];
    for (int c = 0; c < n; ++c)
    {
      B.setEntry(p, c, _b[v_p][c]);
    }
  }
  
  return B;
}

int PhylogeneticTree::getParallelEvolutionCount() const
{
  const int n = _b[_root].size();
  StlIntVector parallelEvolution(n, 0);
  
  for (ArcIt a(_T); a != lemon::INVALID; ++a)
  {
    for (const IntPair& ci : _charStateLabeling[a])
    {
      if (ci.second == 1)
      {
        ++parallelEvolution[ci.first];
      }
    }
  }
  
  int res = 0;
  for (int c = 0; c < n; ++c)
  {
    if (parallelEvolution[c] > 1)
      ++res;
  }
  
  return res;
}

int PhylogeneticTree::getBackMutationCount() const
{
  const int n = _b[_root].size();
  StlIntVector backMutation(n, 0);
  
  for (ArcIt a(_T); a != lemon::INVALID; ++a)
  {
    for (const IntPair& ci : _charStateLabeling[a])
    {
      if (ci.second == 0 || ci.second >= 2)
      {
        ++backMutation[ci.first];
      }
    }
  }
  
  int res = 0;
  for (int c = 0; c < n; ++c)
  {
    if (backMutation[c] > 0)
      ++res;
  }
  
  return res;
}

std::ostream& operator<<(std::ostream& out, const PhylogeneticTree& T)
{
  const int nrNodes = lemon::countNodes(T._T);
  out << nrNodes << " #nodes" << std::endl;
  
  IntNodeMap nodeToIndex(T._T, -1);
  NodeVector indexToNode;
  for (NodeIt v(T._T); v != lemon::INVALID; ++v)
  {
    nodeToIndex[v] = indexToNode.size();
    indexToNode.push_back(v);
  }
  
  bool first = true;
  for (NodeIt v(T._T); v != lemon::INVALID; ++v)
  {
    if (first)
      first = false;
    else
      out << " ";
    
    if (v == T._root)
    {
      out << "-1";
    }
    else
    {
      Node u = T._T.source(InArcIt(T._T, v));
      out << nodeToIndex[u];
    }
  }
  out << std::endl;
  
  for (NodeIt v(T._T); v != lemon::INVALID; ++v)
  {
    if (OutArcIt(T._T, v) == lemon::INVALID)
    {
      out << T._leafToTaxon[v];
    }
    else
    {
      out << "v" << T._T.id(v);
    }
    
    const StlIntVector& b_v = T._b[v];
    const int n = b_v.size();
    for (int c = 0; c < n; ++c)
    {
      out << " " << b_v[c];
    }
    out << std::endl;
  }
  
  return out;
}

std::istream& operator>>(std::istream& in, PhylogeneticTree& T)
{
  g_lineNumber = 0;
  
  std::string line;
  getline(in, line);
  
  std::stringstream ss(line);
  
  int nrNodes = -1;
  ss >> nrNodes;
  
  if (nrNodes < 0)
  {
    throw std::runtime_error(getLineNumber()
                             + "Error: number of nodes should be positive.");
  }
  
  IntNodeMap nodeToIndex(T._T, -1);
  NodeVector indexToNode(nrNodes, lemon::INVALID);
  T._T.clear();
  for (int nodeIndex = 0; nodeIndex < nrNodes; ++nodeIndex)
  {
    Node v = T._T.addNode();
    nodeToIndex[v] = nodeIndex;
    indexToNode[nodeIndex] = v;
  }
  
  getline(in, line);
  ss.clear();
  ss.str(line);

  StlIntVector pi(nrNodes, -2);
  for (int nodeIndex = 0; nodeIndex < nrNodes; ++nodeIndex)
  {
    ss >> pi[nodeIndex];
    
    Node child = indexToNode[nodeIndex];
    if (pi[nodeIndex] != -1)
    {
      Node parent = indexToNode[pi[nodeIndex]];
      T._T.addArc(parent, child);
    }
    else
    {
      T._root = indexToNode[nodeIndex];
    }
  }
  
  // identify leaves
  int nrLeaves = 0;
  for (int nodeIndex = 0; nodeIndex < nrNodes; ++nodeIndex)
  {
    Node v = indexToNode[nodeIndex];
    if (OutArcIt(T._T, v) == lemon::INVALID)
    {
      ++nrLeaves;
    }
  }
  
  T._taxonToLeaf = NodeVector(nrLeaves, lemon::INVALID);
  for (int nodeIndex = 0; nodeIndex < nrNodes; ++nodeIndex)
  {
    getline(in, line);
    ss.clear();
    ss.str(line);
    
    Node v = indexToNode[nodeIndex];
    if (OutArcIt(T._T, v) == lemon::INVALID)
    {
      int p = -1;
      ss >> p;
      if (p < 0 || p >= nrLeaves)
      {
        throw std::runtime_error(getLineNumber()
                                 + "Error: invalid taxon index.");
      }
      T._leafToTaxon[v] = p;
      T._taxonToLeaf[p] = v;
    }
    
    std::vector<std::string> s;
    boost::split(s, line, boost::is_any_of(" "));
    T._b[v] = StlIntVector(s.size() - 1, 0);
    for (int idx = 1; idx < s.size(); ++idx)
    {
      int val = boost::lexical_cast<int>(s[idx]);
      if (val != 0 && val != 1)
      {
        throw std::runtime_error(getLineNumber()
                                 + "Error: invalid state.");
      }
      T._b[v][idx-1] = val;
    }
  }
  
  return in;
}
