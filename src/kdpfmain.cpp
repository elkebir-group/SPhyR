/*
 * kdpfmain.cpp
 *
 *  Created on: 8-mar-2018
 *      Author: M. El-Kebir
 */

#include <fstream>
#include <lemon/arg_parser.h>
#include "matrix.h"
#include "ilpsolverdolloflip.h"
#include "phylogenetictree.h"
#include "columngenflip.h"

int main(int argc, char** argv)
{
  std::string filename;
  int k = 1;
  int memoryLimit = -1;
  int nrThreads = 1;
  int timeLimit = -1;
  double alpha = 1e-3;
  double beta = 0.3;
  bool verbose = false;
  bool columnGeneration = false;
  bool lazy = false;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("c", "Enable column generation", columnGeneration)
    .refOption("k", "Maximum number of losses per SNV (default: 1)", k)
    .refOption("T", "Time limit in seconds (default: -1, unlimited)", timeLimit)
    .refOption("t", "Number of threads (default: 1)", nrThreads)
    .refOption("M", "Memory limit in MB (default: -1, unlimited)", memoryLimit)
    .refOption("v", "Verbose output", verbose)
    .refOption("a", "False positive rate (default: 1e-3)", alpha)
    .refOption("b", "False negative rate (default: 0.3)", beta)
    .refOption("lazy", "Use lazy constraints", lazy)
    .other("input", "Input file")
    .other("output", "Output file");
  ap.parse();
  
  if (ap.files().empty())
  {
    std::cerr << "Error: missing input file" << std::endl;
    return 1;
  }
  
  std::ifstream inD(ap.files()[0]);
  if (!inD.good())
  {
    std::cerr << "Error: failed to open '" << argv[1] << "' for reading"
    << std::endl;
    return 1;
  }
  
  std::string outputFilename = ap.files().size() > 1 ? ap.files()[1] : "";
  
  Matrix D;
  inD >> D;
  inD.close();
  
  StlIntVector mapping;
  D = D.simplify(mapping);
  
  if (columnGeneration)
  {
    ColumnGenFlip solver(D, k, lazy, alpha, beta);
    solver.init();
    
    if (solver.solve(timeLimit, memoryLimit, nrThreads, verbose))
    {
      Matrix A = solver.getSolA().expand(mapping);
      if (outputFilename.empty())
      {
        std::cout << A;
      }
      else
      {
        std::ofstream outE(outputFilename.c_str());
        outE << A;
        outE.close();
      }
    }
  }
  else
  {
    IlpSolverDolloFlip solver(D, k, alpha, beta);
    solver.init();
    
    if (solver.solve(timeLimit, memoryLimit, nrThreads, verbose))
    {
      Matrix A = solver.getSolE().expand(mapping);
      if (outputFilename.empty())
      {
        std::cout << A;
      }
      else
      {
        std::ofstream outE(outputFilename.c_str());
        outE << A;
        outE.close();
      }
    }
  }
  
  return 0;
}
