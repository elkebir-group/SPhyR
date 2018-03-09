/*
 * callback.cpp
 *
 *  Created on: 23-feb-2018
 *      Author: M. El-Kebir
 */

#include "callback.h"

void UserCutCallback::main()
{
  IloNumArray vals(getEnv(), _vars.getSize());
  getValues(vals, _vars);
  
//  int _thresh = std::numeric_limits<int>::max();
  int _thresh = 100;
  
//  std::cout << "Hello" << std::endl;
  
  // TODO this is probably to slow. Can we save the lists?
  int counter = 0;
  for (int j = 0; j < _n && counter < _thresh; j++)
  {
    for (int k = j + 1; k < _n && counter < _thresh; k++)
    {
      // TODO: generalize
      if (j / 2 == k / 2) continue;
      
      std::vector<std::pair<double, int> > list1, list2, list3;
      for (int i = 0; i < _m; i++)
      {
        list1.push_back(std::make_pair(vals[i*_n + j] - vals[i*_n + k], i));
        list2.push_back(std::make_pair(vals[i*_n + k] - vals[i*_n + j], i));
        list3.push_back(std::make_pair(vals[i*_n + j] + vals[i*_n + k], i));
      }
      
      std::sort(list1.begin(), list1.end());
      std::sort(list2.begin(), list2.end());
      std::sort(list3.begin(), list3.end());
      
      int t1, t2, t3;
      double x1, x2, x3;
      for (int i1 = list1.size() - 1; i1 >= 0 && counter < _thresh; i1--)
      {
        t1 = list1[i1].second;
        x1 = list1[i1].first;
        for (int i2 = list2.size() - 1; i2 >= 0 && counter < _thresh; i2--)
        {
          t2 = list2[i2].second;
          x2 = list2[i2].first;
          for (int i3 = list3.size() - 1; i3 >= 0 && counter < _thresh; i3--)
          {
            t3 = list3[i3].second;
            x3 = list3[i3].first;
            if (g_tol.less(3., x1 + x2 + x3))
            {
              counter++;
              add(_vars[t1*_n + j] + _vars[t2*_n + k] + _vars[t3*_n + j] + _vars[t3*_n + k]
                  - _vars[t1*_n + k] - _vars[t2*_n + j] <= 3, IloCplex::UseCutPurge).end();
            }
          }
        }
      }
    }
  }
  
//  std::cout << "Separated (user): " << counter << std::endl;
}

void LazyCutCallback::main()
{
  IloNumArray vals(getEnv(), _vars.getSize());
  getValues(vals, _vars);
  
  int _thresh = std::numeric_limits<int>::max();
//  int _thresh = 100;

  //  std::cout << "Hello" << std::endl;
  
  // TODO this is probably to slow. Can we save the lists?
  int counter = 0;
  for (int j = 0; j < _n && counter < _thresh; j++)
  {
    for (int k = j + 1; k < _n && counter < _thresh; k++)
    {
      // TODO: generalize
      if (j / 2 == k / 2) continue;

      std::vector<std::pair<double, int> > list1, list2, list3;
      for (int i = 0; i < _m; i++)
      {
        list1.push_back(std::make_pair(vals[i*_n + j] - vals[i*_n + k], i));
        list2.push_back(std::make_pair(vals[i*_n + k] - vals[i*_n + j], i));
        list3.push_back(std::make_pair(vals[i*_n + j] + vals[i*_n + k], i));
      }
      
      std::sort(list1.begin(), list1.end());
      std::sort(list2.begin(), list2.end());
      std::sort(list3.begin(), list3.end());
      
      int t1, t2, t3;
      double x1, x2, x3;
      for (int i1 = list1.size() - 1; i1 >= 0 && counter < _thresh; i1--)
      {
        t1 = list1[i1].second;
        x1 = list1[i1].first;
        for (int i2 = list2.size() - 1; i2 >= 0 && counter < _thresh; i2--)
        {
          t2 = list2[i2].second;
          x2 = list2[i2].first;
          for (int i3 = list3.size() - 1; i3 >= 0 && counter < _thresh; i3--)
          {
            t3 = list3[i3].second;
            x3 = list3[i3].first;
            if (g_tol.less(3., x1 + x2 + x3))
            {
              counter++;
              add(_vars[t1*_n + j] + _vars[t2*_n + k] + _vars[t3*_n + j] + _vars[t3*_n + k]
                  - _vars[t1*_n + k] - _vars[t2*_n + j] <= 3, IloCplex::UseCutPurge).end();
            }
          }
        }
      }
    }
  }
  
//  std::cout << "Separated (lazy): " << counter << std::endl;
}
