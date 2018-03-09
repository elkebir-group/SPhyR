/*
 * dollocallback.h
 *
 *  Created on: 25-feb-2018
 *      Author: M. El-Kebir
 */

#ifndef DOLLOCALLBACK_H
#define DOLLOCALLBACK_H

#include <ilcplex/ilocplex.h>
#include "utils.h"

template<class T>
class DolloCallback : public T
{
private:
  const int _m;
  const int _n;
  const int _k;
  IloBoolVarArray _vars;
  IloNumArray _vals;
  int _maxIterations;
  int _currentIterations;
  IloCplex::MIPCallbackI::NodeId _nodeId;
  
  int getIndex(int p, int c, int i) const
  {
    return (_n * (_k + 2)) * p + (_k + 2) * c + i;
  }
  
public:
  typedef IloArray<IloArray<IloBoolVarArray> > IloBoolVar3Matrix;
  
  DolloCallback(IloEnv env,
                const IloBoolVar3Matrix& E,
                const int m,
                const int n,
                const int k)
    : T(env)
    , _m(m)
    , _n(n)
    , _k(k)
    , _vars()
    , _vals()
    , _maxIterations(100)
    , _currentIterations(0)
    , _nodeId()
  {
    _vars = IloBoolVarArray(env, _m * _n * (_k + 2));
    
    for (int p = 0; p < _m; ++p)
    {
      for (int c = 0; c < _n; ++c)
      {
        for (int i = 0; i <= _k + 1; ++i)
        {
          _vars[getIndex(p, c, i)] = E[p][c][i];
        }
      }
    }
    
    _vals = IloNumArray(env, _vars.getSize());
  }
  
  IloCplex::CallbackI *duplicateCallback() const
  {
    return (new (T::getEnv()) DolloCallback(*this));
  }
  
  void main();
  
  void separate();
};

template<class T>
void DolloCallback<T>::separate()
{
  T::getValues(_vals, _vars);
  
  for (int c = 0; c < _n; ++c)
  {
    for (int d = c + 1; d < _n; ++d)
    {
      // i == 1 && j == 1
      int p_star = -1;
      double val_p_star = 0;
      for (int p = 0; p < _m; p++)
      {
        double val = (1 - _vals[getIndex(p, c, 0)]) + _vals[getIndex(p, d, 0)];
        if (val > val_p_star)
        {
          val_p_star = val;
          p_star = p;
        }
      }
      
      int q_star = -1;
      double val_q_star = 0;
      for (int q = 0; q < _m; q++)
      {
        if (q == p_star) continue;
        double val = _vals[getIndex(q, c, 0)] + (1 - _vals[getIndex(q, d, 0)]);
        if (val > val_q_star)
        {
          val_q_star = val;
          q_star = q;
        }
      }
      
      int r_star = -1;
      double val_r_star = 0;
      for (int r = 0; r < _m; r++)
      {
        if (r == p_star) continue;
        if (r == q_star) continue;
        double val = (1 - _vals[getIndex(r, c, 0)]) + (1 - _vals[getIndex(r, d, 0)]);
        if (val > val_r_star)
        {
          val_r_star = val;
          r_star = r;
        }
      }
      
      if (g_tol.less(5., val_p_star + val_q_star + val_r_star))
      {
        assert(p_star != q_star);
        assert(p_star != r_star);
        assert(q_star != r_star);
        assert(g_tol.less(5., (1 - _vals[getIndex(p_star, c, 0)]) + _vals[getIndex(p_star, d, 0)] +
                          _vals[getIndex(q_star, c, 0)] + (1 - _vals[getIndex(q_star, d, 0)]) +
                          (1 - _vals[getIndex(r_star, c, 0)]) + (1 - _vals[getIndex(r_star, d, 0)])));
        
        T::add((1 - _vars[getIndex(p_star, c, 0)]) + _vars[getIndex(p_star, d, 0)] +
               _vars[getIndex(q_star, c, 0)] + (1 - _vars[getIndex(q_star, d, 0)]) +
               (1 - _vars[getIndex(r_star, c, 0)]) + (1 - _vars[getIndex(r_star, d, 0)]) <= 5,
               IloCplex::UseCutPurge).end();
      }
      
      // i == 1, 2 <= j <= _k + 1
      for (int j = 2; j <= _k + 1; ++j)
      {
        p_star = -1;
        val_p_star = 0;
        for (int p = 0; p < _m; p++)
        {
          double val = (1 - _vals[getIndex(p, c, 0)]) + (1 - _vals[getIndex(p, d, j)]);
          if (val > val_p_star)
          {
            val_p_star = val;
            p_star = p;
          }
        }
        
        q_star = -1;
        val_q_star = 0;
        for (int q = 0; q < _m; q++)
        {
          if (q == p_star) continue;
          double val = _vals[getIndex(q, c, 0)] + _vals[getIndex(q, d, j)];
          if (val > val_q_star)
          {
            val_q_star = val;
            q_star = q;
          }
        }
        
        r_star = -1;
        val_r_star = 0;
        for (int r = 0; r < _m; r++)
        {
          if (r == p_star) continue;
          if (r == q_star) continue;
          double val = (1 - _vals[getIndex(r, c, 0)]) + _vals[getIndex(r, d, j)];
          if (val > val_r_star)
          {
            val_r_star = val;
            r_star = r;
          }
        }
        
        if (g_tol.less(5., val_p_star + val_q_star + val_r_star))
        {
          assert(p_star != q_star);
          assert(p_star != r_star);
          assert(q_star != r_star);
          
          assert(g_tol.less(5., (1 - _vals[getIndex(p_star, c, 0)]) + (1 - _vals[getIndex(p_star, d, j)]) +
                            _vals[getIndex(q_star, c, 0)] + _vals[getIndex(q_star, d, j)] +
                            (1 - _vals[getIndex(r_star, c, 0)]) + _vals[getIndex(r_star, d, j)]));
          
          T::add((1 - _vars[getIndex(p_star, c, 0)]) + (1 - _vars[getIndex(p_star, d, j)]) +
                 _vars[getIndex(q_star, c, 0)] + _vars[getIndex(q_star, d, j)] +
                 (1 - _vars[getIndex(r_star, c, 0)]) + _vars[getIndex(r_star, d, j)] <= 5,
                 IloCplex::UseCutPurge).end();
        }
      }
      
      // 2 <= i <= _k + 1, j = 1
      for (int i = 2; i <= _k + 1; ++i)
      {
        p_star = -1;
        val_p_star = 0;
        for (int p = 0; p < _m; p++)
        {
          double val = _vals[getIndex(p, c, i)] + _vals[getIndex(p, d, 0)];
          if (val > val_p_star)
          {
            val_p_star = val;
            p_star = p;
          }
        }
        
        q_star = -1;
        val_q_star = 0;
        for (int q = 0; q < _m; q++)
        {
          if (q == p_star) continue;
          double val = (1 - _vals[getIndex(q, c, i)]) + (1 - _vals[getIndex(q, d, 0)]);
          if (val > val_q_star)
          {
            val_q_star = val;
            q_star = q;
          }
        }
        
        r_star = -1;
        val_r_star = 0;
        for (int r = 0; r < _m; r++)
        {
          if (r == p_star) continue;
          if (r == q_star) continue;
          double val = _vals[getIndex(r, c, i)] + (1 - _vals[getIndex(r, d, 0)]);
          if (val > val_r_star)
          {
            val_r_star = val;
            r_star = r;
          }
        }
        
        if (g_tol.less(5., val_p_star + val_q_star + val_r_star))
        {
          assert(p_star != q_star);
          assert(p_star != r_star);
          assert(q_star != r_star);
          
          assert(g_tol.less(5., _vals[getIndex(p_star, c, i)] + _vals[getIndex(p_star, d, 0)] +
                            (1 - _vals[getIndex(q_star, c, i)]) + (1 - _vals[getIndex(q_star, d, 0)]) +
                            _vals[getIndex(r_star, c, i)] + (1 - _vals[getIndex(r_star, d, 0)])));
          
          T::add(_vars[getIndex(p_star, c, i)] + _vars[getIndex(p_star, d, 0)] +
                 (1 - _vars[getIndex(q_star, c, i)]) + (1 - _vars[getIndex(q_star, d, 0)]) +
                 _vars[getIndex(r_star, c, i)] + (1 - _vars[getIndex(r_star, d, 0)]) <= 5,
                 IloCplex::UseCutPurge).end();
        }
      }
      
      // 2 <= i <= _k + 1, 2 <= j <= _k + 1
      for (int i = 2; i <= _k + 1; ++i)
      {
        for (int j = 2; j <= _k + 1; ++j)
        {
          p_star = -1;
          val_p_star = 0;
          for (int p = 0; p < _m; p++)
          {
            double val = _vals[getIndex(p, c, i)] + (1 - _vals[getIndex(p, d, j)]);
            if (val > val_p_star)
            {
              val_p_star = val;
              p_star = p;
            }
          }
          
          q_star = -1;
          val_q_star = 0;
          for (int q = 0; q < _m; q++)
          {
            if (q == p_star) continue;
            double val = (1 - _vals[getIndex(q, c, i)]) + (1 - _vals[getIndex(q, d, j)]);
            if (val > val_q_star)
            {
              val_q_star = val;
              q_star = q;
            }
          }
          
          r_star = -1;
          val_r_star = 0;
          for (int r = 0; r < _m; r++)
          {
            if (r == p_star) continue;
            if (r == q_star) continue;
            double val = _vals[getIndex(r, c, i)] + _vals[getIndex(r, d, j)];
            if (val > val_r_star)
            {
              val_r_star = val;
              r_star = r;
            }
          }
          
          if (g_tol.less(5., val_p_star + val_q_star + val_r_star))
          {
            assert(p_star != q_star);
            assert(p_star != r_star);
            assert(q_star != r_star);
            
            assert(g_tol.less(5., _vals[getIndex(p_star, c, i)] + (1 - _vals[getIndex(p_star, d, j)]) +
                              (1 - _vals[getIndex(q_star, c, i)]) + (1 - _vals[getIndex(q_star, d, j)]) +
                              _vals[getIndex(r_star, c, i)] + _vals[getIndex(r_star, d, j)]));
            
            T::add(_vars[getIndex(p_star, c, i)] + (1 - _vars[getIndex(p_star, d, j)]) +
                   (1 - _vars[getIndex(q_star, c, i)]) + (1 - _vars[getIndex(q_star, d, j)]) +
                   _vars[getIndex(r_star, c, i)] + _vars[getIndex(r_star, d, j)] <= 5,
                   IloCplex::UseCutPurge).end();
          }
        }
      }
    }
  }
}

template<>
void DolloCallback<IloCplex::UserCutCallbackI>::main()
{
  if (_currentIterations == _maxIterations)
  {
    abortCutLoop();
  }
  
  IloCplex::MIPCallbackI::NodeId newId = getNodeId();
  if (_currentIterations == 0 || newId != _nodeId)
  {
    _nodeId = newId;
  }
  _currentIterations++;
  
  separate();
}


template<class T>
void DolloCallback<T>::main()
{
//  static typename T::NodeId nodeId;
//  typename T::NodeId nodeZero;
//  if (nodeId == T::getNodeId() && nodeId != nodeZero)
//  {
//    return;
//  }
//  else
//  {
//    nodeId = T::getNodeId();
//  }
  

  
//  for (int p = 0; p < _m; ++p)
//  {
//    for (int c = 0; c < _n; ++c)
//    {
//      std::cout << "( ";
//      for (int i = 0; i <= _k + 1; ++i)
//      {
//        std::cout << _vals[getIndex(p, c, i)] << " ";
//      }
//      std::cout << ") ";
//    }
//    std::cout << std::endl;
//  }
//  std::cout << std::endl;
  
  separate();
}

#endif // DOLLOCALLBACK_H
