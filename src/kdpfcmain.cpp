/*
 * dollomain.cpp
 *
 *  Created on: 9-mar-2018
 *      Author: M. El-Kebir
 */

#include <fstream>
#include <lemon/arg_parser.h>
#include "matrix.h"
//#include "ilpsolverdolloflipcluster.h"
#include "coordinateascent.h"

int main(int argc, char** argv)
{
  std::string filename;
  int k = 1;
  double alpha = 1e-3;
  double beta = 0.3;
  int lC = 10;
  int lT = 10;
  bool exact = false;
  int seed = 0;
  int memoryLimit = -1;
  int nrThreads = 1;
  int timeLimit = -1;
  bool verbose = false;
  bool lazy = true;
  int restarts = 1;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("k", "Maximum number of losses per SNV (default: 1)", k)
    .refOption("a", "False positive rate (default: 1e-3)", alpha)
    .refOption("b", "False negative rate (default: 0.3)", beta)
    .refOption("lC", "Number of SNV clusters (default: 10)", lC)
    .refOption("lT", "Number of cell clusters (default: 10)", lT)
    .refOption("exact", "Exact algorithm", exact)
    .refOption("N", "Number of restarts (default: 1)", restarts)
    .refOption("s", "Random number generator seed (default: 0)", seed)
    .refOption("T", "Time limit in seconds (default: -1, unlimited). Only used with -exact.", timeLimit)
    .refOption("t", "Number of threads (default: 1)", nrThreads)
    .refOption("M", "Memory limit in MB (default: -1, unlimited)", memoryLimit)
    .refOption("v", "Verbose output", verbose)
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
    std::cerr << "Error: failed to open '" << ap.files()[0] << "' for reading"
              << std::endl;
    return 1;
  }
  
  std::string outputFilename = ap.files().size() > 1 ? ap.files()[1] : "";
  
  Matrix D;
  inD >> D;
  inD.close();

  StlIntVector characterMapping, taxonMapping;
  D = D.simplify(characterMapping, taxonMapping);
  
//  if (exact)
//  {
//    IlpSolverDolloFlipCluster solver(D, k, alpha, beta, l);
//    solver.init();
//    if (solver.solve(timeLimit, memoryLimit, nrThreads, verbose))
//    {
//      Matrix A = solver.getSolE().expand(characterMapping, taxonMapping);
//      if (outputFilename.empty())
//      {
//        std::cout << A;
//      }
//      else
//      {
//        std::ofstream outE(outputFilename.c_str());
//        outE << A;
//        outE.close();
//      }
//    }
//  }
//  else
  {
    CoordinateAscent ca(D, k, lazy, alpha, beta, lT, lC, seed);
    ca.solve(-1, memoryLimit, nrThreads, verbose, restarts);
    Matrix bestA = ca.getE();
//    std::cout << bestA << std::endl;
    bestA = bestA.expandColumns(ca.getZC());
//    std::cout << bestA << std::endl;
    bestA = bestA.expandRows(ca.getZT());
//    std::cout << bestA << std::endl;
    bestA = bestA.expand(characterMapping, taxonMapping);
//    std::cout << bestA << std::endl;
    
    std::cerr << "Solution likelihood: " << ca.getLogLikelihood() << std::endl;
    
    if (outputFilename.empty())
    {
      std::cout << bestA;
    }
    else
    {
      std::ofstream outFile(outputFilename.c_str());
      outFile << bestA;
      outFile.close();
    }
  }
}
