/*
 * dollomain.cpp
 *
 *  Created on: 22-feb-2018
 *      Author: M. El-Kebir
 */

#include <fstream>
#include <lemon/arg_parser.h>
#include "matrix.h"
#include "ilpsolverdollo.h"
#include "ilpsolverdolloflip.h"
#include "ilpsolverdolloflipcluster.h"
#include "coordinateascent.h"
#include "phylogenetictree.h"

int main(int argc, char** argv)
{
  std::string filename;
  int k = 1;
  double alpha = 1e-3;
  double beta = 0.1;
  int l = -1;
  bool coordinateAscent = false;
  int seed = 0;
  int memoryLimit = -1;
  int nrThreads = 1;
  int timeLimit = -1;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("k", "Maximum number of losses per SNV (default: 1)", k)
    .refOption("a", "False positive rate (default: 1e-3)", alpha)
    .refOption("b", "False negative rate (default: 0.1)", beta)
    .refOption("l", "Number of SNV clusters (default: -1)", l)
    .refOption("cA", "Enable coordinate ascent", coordinateAscent)
    .refOption("s", "Random number generator seed", seed)
    .refOption("T", "Time limit in seconds (default: -1, unlimited)", timeLimit)
    .refOption("t", "Number of threads (default: 1)", nrThreads)
    .refOption("M", "Memory limit in MB (default: -1, unlimited)", memoryLimit)
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
  
  std::cout << D << std::endl;
  
  IlpSolverDollo* pSolver = NULL;
  if (g_tol.nonZero(alpha) || g_tol.nonZero(beta))
  {
    if (l == -1)
    {
      pSolver = new IlpSolverDolloFlip(D, k, alpha, beta);
    }
    else
    {
      if (!coordinateAscent)
      {
        pSolver = new IlpSolverDolloFlipCluster(D, k, alpha, beta, l);
      }
      else
      {
        CoordinateAscent ca(D, k, alpha, beta, l, seed);
        ca.solve(timeLimit, memoryLimit, nrThreads);
        
        if (outputFilename.empty())
        {
          PhylogeneticTree T(ca.getE(), k);
          if (T.reconstructTree())
          {
            T.writeDOT(std::cout);
          }
        }
        else
        {
          std::ofstream outFile(outputFilename.c_str());
          outFile << ca.getE();
          outFile << ca.getZ();
          outFile.close();
        }
        
        return 0;
      }
    }
  }
  else
  {
    pSolver = new IlpSolverDollo(D, k);
  }
  
  pSolver->init();
  pSolver->solve(timeLimit, memoryLimit, nrThreads);
  delete pSolver;
  
  return 0;
}
