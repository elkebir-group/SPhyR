/*
 * simulatemain.cpp
 *
 *  Created on: 7-mar-2018
 *      Author: M. El-Kebir
 */

#include <lemon/arg_parser.h>
#include "phylogenetictree.h"
#include <fstream>

int main(int argc, char** argv)
{
  int k = -1;
  int seed = 0;
  double lossRate;
  std::string filenameA;
  std::string filenameB;
  std::string filenameDOT;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("k", "Maximum number of losses per SNV", k, true)
    .refOption("s", "Random number generator seeding (default: 0)", seed)
    .refOption("loss", "Loss rate", lossRate, true)
    .refOption("A", "Matrix A output filename", filenameA)
    .refOption("B", "Matrix B output filename", filenameB)
    .refOption("dot", "DOT ouput filename", filenameDOT)
    .other("input", "Input file");
  ap.parse();
  
  if (ap.files().empty())
  {
    std::cerr << "Error: input file missing" << std::endl;
    return 1;
  }
  
  Matrix inputB;
  std::ifstream inB(ap.files()[0]);
  if (!inB.good())
  {
    std::cerr << "Error: could not open '" << ap.files()[0] << "' for reading" << std::endl;
    return 1;
  }
  inB >> inputB;
  inB.close();
  
  PhylogeneticTree phyloT(inputB, 0);
  if (!phyloT.reconstructTree())
  {
    std::cerr << "Error: provided matrix is not a perfect phylogeny matrix" << std::endl;
    return 1;
  }
  
  phyloT.generateLosses(lossRate, k, seed);
  
  if (!filenameDOT.empty())
  {
    std::ofstream outDOT(filenameDOT.c_str());
    phyloT.writeDOT(outDOT);
    outDOT.close();
  }

  if (!filenameA.empty())
  {
    std::ofstream outA(filenameA.c_str());
    outA << phyloT.getA();
    outA.close();
  }

  if (!filenameB.empty())
  {
    std::ofstream outB(filenameB.c_str());
    outB << phyloT.getB();
    outB.close();
  }
  
  if (filenameDOT.empty() && filenameA.empty() && filenameB.empty())
  {
    std::cout << phyloT.getB();
  }
  
  return 0;
}
