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
#include "dollophylogenetictree.h"

int main(int argc, char** argv)
{
  double alpha = 1e-3;
  double beta = 0.3;
  bool tree = false;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("a", "False positive rate (default: 1e-3)", alpha)
    .refOption("b", "False negative rate (default: 0.3)", beta)
    .refOption("T", "Use tree instead of matrix", tree)
    .other("inferred", "Inferred solution file")
    .other("true", "True solution file")
    .other("input", "Input matrix");
  
  ap.parse();
  
  if (ap.files().empty())
  {
    std::cerr << "Error: input file missing" << std::endl;
    return 1;
  }
  
  if (ap.files().size() == 1)
  {
    Matrix::ViolationList list;
    Matrix* pInferredA = Matrix::parse(ap.files()[0]);
    if (!pInferredA)
    {
      return 1;
    }
    pInferredA->identifyViolations(pInferredA->getMaxNrLosses(), list);
    
    IntPairSet violationEntries;
    for (const Matrix::Violation& violation : list)
    {
      std::cout << "Condition " << violation._condition
                << " ; " << "p = " << violation._p << " ; q = " << violation._q
                << " ; r = " << violation._r << " ; c = " << violation._c
                << " ; d = " << violation._d << std::endl;
      std::cout << pInferredA->getEntry(violation._p, violation._c)
                << " " << pInferredA->getEntry(violation._p, violation._d)
                << std::endl;
      std::cout << pInferredA->getEntry(violation._q, violation._c)
                << " " << pInferredA->getEntry(violation._q, violation._d)
                << std::endl;
      std::cout << pInferredA->getEntry(violation._r, violation._c)
                << " " << pInferredA->getEntry(violation._r, violation._d)
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
              << " / " << pInferredA->getNrTaxa() * pInferredA->getNrCharacters()
              << " = " << (double) violationEntries.size() /
                (pInferredA->getNrTaxa() * pInferredA->getNrCharacters())
              << std::endl;
    
    for (int c = 0; c < pInferredA->getNrCharacters(); ++c)
    {
      int ones_c = pInferredA->getNrOfOnes(c);
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
  else if (ap.files().size() == 3)
  {
    PhylogeneticTree* pInferredT = NULL;
    if (tree)
    {
      pInferredT = PhylogeneticTree::parse(ap.files()[0]);
    }
    else
    {
      pInferredT = DolloPhylogeneticTree::parse(ap.files()[0]);
    }
    
    if (!pInferredT)
    {
      return 1;
    }
    
    PhylogeneticTree* pTrueT = DolloPhylogeneticTree::parse(ap.files()[1]);
    if (!pTrueT)
    {
      return 1;
    }
    
    Matrix* pInput = Matrix::parse(ap.files()[2]);
    if (!pInput)
    {
      return 1;
    }
    
    Comparison compare(*pTrueT, *pInferredT);
    
    double ancestralRecall, incomparableRecall, clusteredRecall;
    compare.recallCharStatePairs(ancestralRecall, incomparableRecall, clusteredRecall);
    
    double taxaRI, taxaRecall, taxaPrecision;
    compare.getTaxaClusteringMetrics(taxaRI, taxaRecall, taxaPrecision);
    
    double charactersRI, charactersRecall, charactersPrecision;
    compare.getCharactersClusteringMetrics(charactersRI, charactersRecall, charactersPrecision);
    
    Matrix trueA = pTrueT->getMatrix();
    Matrix inferredA = pInferredT->getMatrix();

    double logLikelihood = pInput->getLogLikelihood(inferredA, alpha, beta);
    
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
    // 15. likelihood
    // 16. back mutation count (inferred)
    // 17. parallel evolution count (inferred)
    // 16. back mutation count (true)
    // 17. parallel evolution count (true)
    std::cout << inferredA.getMaxNrLosses() << ","
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
              << charactersPrecision << ","
              << logLikelihood << ","
              << pInferredT->getBackMutationCount() << ","
              << pInferredT->getParallelEvolutionCount() << ","
              << pTrueT->getBackMutationCount() << ","
              << pTrueT->getParallelEvolutionCount()
              << std::endl;
  }
  
  return 0;
}
