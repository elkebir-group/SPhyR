/*
 * utils.h
 *
 *  Created on: 22-feb-2018
 *      Author: M. El-Kebir
 */

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <lemon/list_graph.h>
#include <lemon/tolerance.h>
#include <lemon/time_measure.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

typedef lemon::ListDigraph Digraph;
DIGRAPH_TYPEDEFS(Digraph);
typedef std::vector<Node> NodeVector;
typedef std::vector<bool> StlBoolVector;
typedef std::vector<double> StlDoubleVector;
typedef std::vector<int> StlIntVector;
typedef std::vector<StlIntVector> StlIntMatrix;
typedef std::set<int> StlIntSet;
typedef std::vector<std::string> StringVector;
typedef std::pair<int, int> IntPair;
typedef std::set<IntPair> IntPairSet;
typedef std::vector<IntPair> IntPairVector;

/// Platform independent extraction of line from input string
///
/// @param is Input stream
/// @param t Output string
std::istream& getline(std::istream& is, std::string& t);

/// Return a string containing the current line number
std::string getLineNumber();

/// Tolerance for floating point comparisons
extern lemon::Tolerance<double> g_tol;

/// Current line number
extern int g_lineNumber;

/// Global timer
extern lemon::Timer g_timer;

#endif // UTILS_H
