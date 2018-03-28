/*
 * visualizemain.cpp
 *
 *  Created on: 28-feb-2018
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "matrix.h"
#include "dollophylogenetictree.h"
#include <lemon/arg_parser.h>
#include <fstream>

StringVector parseLabels(const std::string& filename)
{
  StringVector labels;
  
  std::ifstream inFile(filename.c_str());
  while (inFile.good())
  {
    std::string label;
    getline(inFile, label);
    labels.push_back(label);
  }
  
  return labels;
}

int main(int argc, char** argv)
{
  std::string filename, charLabelsFilename, taxonLabelsFilename;
  bool tree = false;
  lemon::ArgParser ap(argc, argv);
  ap.refOption("c", "Character labels", charLabelsFilename, false)
    .refOption("t", "Taxon labels", taxonLabelsFilename, false)
    .refOption("T", "Use tree instead of matrix", tree)
    .other("input", "Input file");
  ap.parse();
  
  if (ap.files().empty())
  {
    std::cerr << "Error: missing input file" << std::endl;
    return 1;
  }

  StringVector characterLabels = parseLabels(charLabelsFilename);
  StringVector taxonLabels = parseLabels(taxonLabelsFilename);

  if (tree)
  {
    PhylogeneticTree* pTree = PhylogeneticTree::parse(ap.files()[0]);
    if (pTree)
    {
      if (characterLabels.empty() || taxonLabels.empty())
      {
        pTree->writeDOT(std::cout);
      }
      else
      {
        pTree->writeDOT(std::cout, taxonLabels, characterLabels);
      }
    }
    else
    {
      return 1;
    }
  }
  else
  {
    DolloPhylogeneticTree* pTree = DolloPhylogeneticTree::parse(ap.files()[0]);
    if (pTree)
    {
      if (characterLabels.empty() || taxonLabels.empty())
      {
        pTree->writeDOT(std::cout);
      }
      else
      {
        pTree->writeDOT(std::cout, taxonLabels, characterLabels);
      }
    }
    else
    {
      return 1;
    }
  }
    
  return 0;
}
