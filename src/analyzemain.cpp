/*
 * analyzemain.cpp
 *
 *  Created on: 8-mar-2018
 *      Author: M. El-Kebir
 */

#include <lemon/arg_parser.h>
#include <fstream>
#include "matrix.h"
#include "comparison.h"

bool parse(const std::string& filename,
           Matrix& A)
{
  std::ifstream inA(filename.c_str());
  if (!inA.good())
  {
    std::cerr << "Error: could not open '" << filename << "' for reading" << std::endl;
    return false;
  }
  inA >> A;
  inA.close();
  
  return true;
}

int main(int argc, char** argv)
{
  lemon::ArgParser ap(argc, argv);
  ap.other("inferred", "Inferred solution file")
    .other("true", "True solution file");
  
  ap.parse();
  
  if (ap.files().empty())
  {
    std::cerr << "Error: input file missing" << std::endl;
    return 1;
  }
  
  Matrix inferredA;
  if (!parse(ap.files()[0], inferredA))
    return 1;
  
  if (ap.files().size() == 1)
  {
    Matrix::ViolationList list;
    inferredA.identifyViolations(inferredA.getMaxNrLosses(), list);
    
    IntPairSet violationEntries;
    for (const Matrix::Violation& violation : list)
    {
      std::cout << "Condition " << violation._condition
                << " ; " << "p = " << violation._p << " ; q = " << violation._q
                << " ; r = " << violation._r << " ; c = " << violation._c
                << " ; d = " << violation._d << std::endl;
      std::cout << inferredA.getEntry(violation._p, violation._c)
                << " " << inferredA.getEntry(violation._p, violation._d)
                << std::endl;
      std::cout << inferredA.getEntry(violation._q, violation._c)
                << " " << inferredA.getEntry(violation._q, violation._d)
                << std::endl;
      std::cout << inferredA.getEntry(violation._r, violation._c)
                << " " << inferredA.getEntry(violation._r, violation._d)
                << std::endl;
      std::cout << std::endl;
      violationEntries.insert(IntPair(violation._p, violation._c));
      violationEntries.insert(IntPair(violation._q, violation._c));
      violationEntries.insert(IntPair(violation._r, violation._c));
      violationEntries.insert(IntPair(violation._p, violation._d));
      violationEntries.insert(IntPair(violation._q, violation._d));
      violationEntries.insert(IntPair(violation._r, violation._d));
    }
    
    std::cout << "Total number of violations: " << list.size() << std::endl;
    std::cout << "Violated entries: " << violationEntries.size()
              << " / " << inferredA.getNrTaxa() * inferredA.getNrCharacters()
              << " = " << (double) violationEntries.size() /
                (inferredA.getNrTaxa() * inferredA.getNrCharacters())
              << std::endl;
    
    for (int c = 0; c < inferredA.getNrCharacters(); ++c)
    {
      int ones_c = inferredA.getNrOfOnes(c);
      if (ones_c == 0)
      {
        std::cout << "Character " << c << " has no 1s" << std::endl;
      }
      else if (ones_c == 1)
      {
        std::cout << "Character " << c << " has a single 1" << std::endl;
      }
    }
  }
  else
  {
    Matrix trueA;
    if (!parse(ap.files()[1], trueA))
      return 1;
    
    PhylogeneticTree inferredT(inferredA, inferredA.getMaxNrLosses());
    if (!inferredT.reconstructTree())
    {
      std::cerr << "Error: provided matrix '" << ap.files()[0] << "' is not a k-Dollo completion for k = " << inferredA.getMaxNrLosses() << std::endl;
      return 1;
    }
    
    PhylogeneticTree trueT(trueA, trueA.getMaxNrLosses());
    if (!trueT.reconstructTree())
    {
      std::cerr << "Error: provided matrix '" << ap.files()[1] << "' is not a k-Dollo completion for k = " << trueA.getMaxNrLosses() << std::endl;
      return 1;
    }
    
    Comparison compare(trueA, inferredA, trueT, inferredT);
    compare.compare();
    
    double ancestralRecall, incomparableRecall, clusteredRecall;
    compare.recallCharStatePairs(ancestralRecall, incomparableRecall, clusteredRecall);
    
    double taxaRI, taxaRecall, taxaPrecision;
    compare.getTaxaClusteringMetrics(taxaRI, taxaRecall, taxaPrecision);
    
    double charactersRI, charactersRecall, charactersPrecision;
    compare.getCharactersClusteringMetrics(charactersRI, charactersRecall, charactersPrecision);
    
    // 1. filename (inferred)
    // 2. inferredK
    // 3. trueK
    // 4. RF
    // 5. normalized RF
    // 6. ancestral pairs recall
    // 7. incomparable pairs recall
    // 8. clustered pairs recall
    // 9. taxa Rand index
    // 10. taxa recall
    // 11. taxa precision
    // 12. character Rand index
    // 13. character recall
    // 14. character precision
    std::cout << ap.files()[0] << ","
              << inferredA.getMaxNrLosses() << ","
              << trueA.getMaxNrLosses() << ","
              << compare.getRF() << ","
              << compare.getNormalizedRF() << ","
              << ancestralRecall << ","
              << incomparableRecall << ","
              << clusteredRecall << ","
              << taxaRI << ","
              << taxaRecall << ","
              << taxaPrecision << ","
              << charactersRI << ","
              << charactersRecall << ","
              << charactersPrecision
              << std::endl;
  }
  
  
  
  return 0;
}
