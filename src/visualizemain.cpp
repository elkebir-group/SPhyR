/*
 * visualizemain.cpp
 *
 *  Created on: 28-feb-2018
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "matrix.h"
#include "phylogenetictree.h"
#include <lemon/arg_parser.h>
#include <fstream>

int main(int argc, char** argv)
{
  std::string filename, charLabelsFilename, taxonLabelsFilename;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("c", "Character labels", charLabelsFilename, false)
    .refOption("t", "Taxon labels", taxonLabelsFilename, false)
    .other("input", "Input file");
  ap.parse();
  
  if (ap.files().empty())
  {
    std::cerr << "Error: missing input file" << std::endl;
    return 1;
  }
  
  std::ifstream inE(ap.files()[0]);
  if (!inE.good())
  {
    std::cerr << "Error: failed to open '" << ap.files()[0] << "' for reading"
              << std::endl;
    return 1;
  }
  
  Matrix E;
  inE >> E;
  inE.close();
  
  int k = E.getMaxNrLosses();
  PhylogeneticTree T(E, k);
  if (!T.reconstructTree())
  {
    std::cerr << "Error: provided matrix is not a k-Dollo completion for k = " << k << std::endl;
    return 1;
  }
  T.writeDOT(std::cout);
  
  return 0;
}
