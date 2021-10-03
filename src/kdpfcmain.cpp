/*
 * kdpfcmain.cpp
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
  int t = 15;
  int s = 10;
  int seed = 0;
  int memoryLimit = -1;
  int nrThreads = 1;
  int timeLimit = -1;
  bool verbose = false;
  bool lazy = true;
  int restarts = 10;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("k", "Maximum number of losses per SNV (default: 1)", k)
    .refOption("a", "False positive rate (default: 1e-3)", alpha)
    .refOption("b", "False negative rate (default: 0.3)", beta)
    .refOption("t", "Number of character clusters (default: 15)", t)
    .synonym("lC", "t")
    .refOption("s", "Number of taxon clusters (default: 10)", s)
    .synonym("lT", "s")
    .refOption("N", "Number of restarts (default: 10)", restarts)
    .refOption("s", "Random number generator seed (default: 0)", seed)
    .refOption("T", "Time limit in seconds (default: -1, unlimited).", timeLimit)
    .refOption("t", "Number of threads (default: 1)", nrThreads)
    .refOption("M", "Memory limit in MB (default: -1, unlimited)", memoryLimit)
    .refOption("v", "Verbose output", verbose)
    .other("input", "Input file")
    .other("output", "Output file");
  ap.parse();

  Matrix D;
  if (!ap.files().empty())
  {
    std::ifstream inD(ap.files()[0]);
    if (!inD.good())
    {
      std::cerr << "Error: failed to open '" << ap.files()[0] << "' for reading"
                << std::endl;
      return 1;
    }
    
    inD >> D;
    inD.close();
  }
  else
  {
    std::cin >> D;
  }
  
  std::string outputFilename = ap.files().size() > 1 ? ap.files()[1] : "";
    
  StlIntVector characterMapping, taxonMapping;
  Matrix simpleD = D.simplify(characterMapping, taxonMapping);
  
  CoordinateAscent ca(simpleD,
                      characterMapping,
                      taxonMapping,
                      k, lazy, alpha, beta, s, t, seed);
  ca.solve(timeLimit, memoryLimit, nrThreads, verbose, restarts);
  Matrix bestA = ca.getE();
  bestA = bestA.expandColumns(ca.getZC());
  bestA = bestA.expandRows(ca.getZT());
  bestA = bestA.expand(characterMapping, taxonMapping);
  
  std::cerr << "Solution likelihood: " << ca.getLogLikelihood() << std::endl;
  
//  assert(!g_tol.different(ca.getLogLikelihood(), D.getLogLikelihood(bestA, alpha, beta)));
  
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
