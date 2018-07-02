/*
 * kdpmain.cpp
 *
 *  Created on: 8-mar-2018
 *      Author: M. El-Kebir
 */

#include <fstream>
#include <lemon/arg_parser.h>
#include "matrix.h"
#include "phylogenetictree.h"
#include "columngen.h"

int main(int argc, char** argv)
{
  std::string filename;
  int k = 1;
  int memoryLimit = -1;
  int nrThreads = 1;
  int timeLimit = -1;
  bool verbose = false;
  bool lazy = true;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("k", "Maximum number of losses per character (default: 1)", k)
    .refOption("T", "Time limit in seconds (default: -1, unlimited)", timeLimit)
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
  
  StlIntVector chacterMapping, taxonMapping;
  D = D.simplify(chacterMapping, taxonMapping);
  
  ColumnGen solver(D, k, lazy);
  solver.init();
  if (solver.solve(timeLimit, memoryLimit, nrThreads, verbose))
  {
    if (outputFilename.empty())
    {
      std::cout << solver.getSolA().expand(chacterMapping, taxonMapping);
    }
    else
    {
      std::ofstream outE(outputFilename.c_str());
      outE << solver.getSolA().expand(chacterMapping, taxonMapping);
      outE.close();
    }
  }
  else
  {
    return 1;
  }
  
  return 0;
}
