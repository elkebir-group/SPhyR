/*
 * dollocallback.h
 *
 *  Created on: 25-feb-2018
 *      Author: M. El-Kebir
 */

#ifndef DOLLOCALLBACK_H
#define DOLLOCALLBACK_H

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilothread.h>
#include "utils.h"

template<class T>
class DolloCallback : public T
{
private:
  const int _m;
  const int _n;
  const int _k;
  IloBoolVarArray _vars;
  int _maxIterations;
  int _currentIterations;
  IloCplex::MIPCallbackI::NodeId _nodeId;
  IloFastMutex* _pMutex;
  
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
                const int k,
                IloFastMutex* pMutex)
    : T(env)
    , _m(m)
    , _n(n)
    , _k(k)
    , _vars()
    , _maxIterations(100)
    , _currentIterations(0)
    , _nodeId()
    , _pMutex(pMutex)
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
  IloNumArray vals = IloNumArray(T::getEnv(), _vars.getSize());

  _pMutex->lock();
  T::getValues(vals, _vars);
  _pMutex->unlock();
  
  for (int c = 0; c < _n; ++c)
  {
    for (int d = c + 1; d < _n; ++d)
    {
      // condition 1
      for (int j = 2; j <= _k + 1; ++j)
      {
        for (int j_prime = 1; j_prime <= _k + 1; ++j_prime)
        {
          if (j_prime == j) continue;
          
          int p_star = -1;
          double val_p_star = 0;
          for (int p = 0; p < _m; p++)
          {
            double val = vals[getIndex(p, c, 1)] + vals[getIndex(p, d, j_prime)];
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
            double val = vals[getIndex(q, c, 0)] + vals[getIndex(q, d, j)];
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
            double val = vals[getIndex(r, c, 1)] + vals[getIndex(r, d, j)];
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
            
            assert(g_tol.less(5., vals[getIndex(p_star, c, 1)] + vals[getIndex(p_star, d, j_prime)] +
                              vals[getIndex(q_star, c, 0)] + vals[getIndex(q_star, d, j)] +
                              vals[getIndex(r_star, c, 1)] + vals[getIndex(r_star, d, j)]));
            
            T::add(_vars[getIndex(p_star, c, 1)] + _vars[getIndex(p_star, d, j_prime)] +
                   _vars[getIndex(q_star, c, 0)] + _vars[getIndex(q_star, d, j)] +
                   _vars[getIndex(r_star, c, 1)] + _vars[getIndex(r_star, d, j)] <= 5,
                   IloCplex::UseCutPurge).end();
          }
        }
      }
      
      // condition 2
      for (int i = 2; i <= _k + 1; ++i)
      {
        for (int i_prime = 1; i_prime <= _k + 1; ++i_prime)
        {
          if (i == i_prime) continue;
          
          int p_star = -1;
          double val_p_star = 0;
          for (int p = 0; p < _m; p++)
          {
            double val = vals[getIndex(p, c, i)] + vals[getIndex(p, d, 0)];
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
            double val = vals[getIndex(q, c, i_prime)] + vals[getIndex(q, d, 1)];
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
            double val = vals[getIndex(r, c, i)] + vals[getIndex(r, d, 1)];
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
            
            assert(g_tol.less(5., vals[getIndex(p_star, c, i)] + vals[getIndex(p_star, d, 0)] +
                              vals[getIndex(q_star, c, i_prime)] + vals[getIndex(q_star, d, 1)] +
                              vals[getIndex(r_star, c, i)] + vals[getIndex(r_star, d, 1)]));
            
            T::add(_vars[getIndex(p_star, c, i)] + _vars[getIndex(p_star, d, 0)] +
                   _vars[getIndex(q_star, c, i_prime)] + _vars[getIndex(q_star, d, 1)] +
                   _vars[getIndex(r_star, c, i)] + _vars[getIndex(r_star, d, 1)] <= 5,
                   IloCplex::UseCutPurge).end();
          }
        }
      }
      
      // condition 3
      for (int i = 2; i <= _k + 1; ++i)
      {
        for (int i_prime = 1; i_prime <= _k + 1; ++i_prime)
        {
          if (i == i_prime) continue;
          for (int j = 2; j <= _k + 1; ++j)
          {
            for (int j_prime = 1; j_prime <= _k + 1; ++j_prime)
            {
              if (j == j_prime) continue;
              
              int p_star = -1;
              double val_p_star = 0;
              for (int p = 0; p < _m; p++)
              {
                double val = vals[getIndex(p, c, i)] + vals[getIndex(p, d, j_prime)];
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
                double val = vals[getIndex(q, c, i_prime)] + vals[getIndex(q, d, j)];
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
                double val = vals[getIndex(r, c, i)] + vals[getIndex(r, d, j)];
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
                
                assert(g_tol.less(5., vals[getIndex(p_star, c, i)] + vals[getIndex(p_star, d, j_prime)] +
                                  vals[getIndex(q_star, c, i_prime)] + vals[getIndex(q_star, d, j)] +
                                  vals[getIndex(r_star, c, i)] + vals[getIndex(r_star, d, j)]));
                
                T::add(_vars[getIndex(p_star, c, i)] + _vars[getIndex(p_star, d, j_prime)] +
                       _vars[getIndex(q_star, c, i_prime)] + _vars[getIndex(q_star, d, j)] +
                       _vars[getIndex(r_star, c, i)] + _vars[getIndex(r_star, d, j)] <= 5,
                       IloCplex::UseCutPurge).end();
              }
            }
          }
        }
      }
      
      // condition 4
      for (int i = 1; i <= _k + 1; ++i)
      {
        for (int i_prime = 1; i_prime <= _k + 1; ++i_prime)
        {
          for (int j = 1; j <= _k + 1; ++j)
          {
            for (int j_prime = 1; j_prime <= _k + 1; ++j_prime)
            {
              int p_star = -1;
              double val_p_star = 0;
              for (int p = 0; p < _m; p++)
              {
                double val = vals[getIndex(p, c, i)] + vals[getIndex(p, d, 0)];
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
                double val = vals[getIndex(q, c, 0)] + vals[getIndex(q, d, j)];
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
                double val = vals[getIndex(r, c, i_prime)] + vals[getIndex(r, d, j_prime)];
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
                assert(g_tol.less(5., vals[getIndex(p_star, c, i)] + vals[getIndex(p_star, d, 0)] +
                                  vals[getIndex(q_star, c, 0)] + vals[getIndex(q_star, d, j)] +
                                  vals[getIndex(r_star, c, i_prime)] + vals[getIndex(r_star, d, j_prime)]));
                
                T::add(_vars[getIndex(p_star, c, i)] + _vars[getIndex(p_star, d, 0)] +
                       _vars[getIndex(q_star, c, 0)] + _vars[getIndex(q_star, d, j)] +
                       _vars[getIndex(r_star, c, i_prime)] + _vars[getIndex(r_star, d, j_prime)] <= 5,
                       IloCplex::UseCutPurge).end();
              }
            }
          }
        }
      }
    }
  }
}

template<>
inline void DolloCallback<IloCplex::UserCutCallbackI>::main()
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
inline void DolloCallback<T>::main()
{
  if (_currentIterations == _maxIterations)
  {
    return;
  }
  
  IloCplex::MIPCallbackI::NodeId newId = T::getNodeId();
  if (_currentIterations == 0 || newId != _nodeId)
  {
    _nodeId = newId;
  }
  _currentIterations++;
  
  separate();
}

#endif // DOLLOCALLBACK_H
