/*
 * perturbmain.cpp
 *
 *  Created on: 8-mar-2018
 *      Author: M. El-Kebir
 */

#include <lemon/arg_parser.h>
#include <fstream>
#include "matrix.h"

int main(int argc, char** argv)
{
  double alpha = 1e-3;
  double beta = 0.3;
  int seed = 0;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("a", "False positive rate (default: 1e-3)", alpha)
    .refOption("b", "False negative rate (default: 0.3)", beta)
    .refOption("s", "Random number generator seed (default: 0)", seed)
    .other("input", "Input file")
    .other("output", "Output file");
  ap.parse();
  
  if (ap.files().empty())
  {
    std::cerr << "Error: input file missing" << std::endl;
    return 1;
  }
  
  Matrix B;
  std::ifstream inB(ap.files()[0]);
  if (!inB.good())
  {
    std::cerr << "Error: could not open '" << ap.files()[0] << "' for reading" << std::endl;
    return 1;
  }
  inB >> B;
  inB.close();
  
  std::string outputFilename = ap.files().size() > 1 ? ap.files()[1] : "";
  
  B.perturb(alpha, beta, seed);
  if (outputFilename.empty())
  {
    std::cout << B;
  }
  else
  {
    std::ofstream outFile(outputFilename.c_str());
    outFile << B;
    outFile.close();
  }
  
  return 0;
}
