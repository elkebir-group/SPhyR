/*
 * analyzemain.cpp
 *
 *  Created on: 8-mar-2018
 *      Author: M. El-Kebir
 */

#include <lemon/arg_parser.h>
#include <fstream>
#include "matrix.h"

int main(int argc, char** argv)
{
  int k = -1;
  
  lemon::ArgParser ap(argc, argv);
  ap.refOption("k", "Maximum number of losses per SNV", k, true)
    .other("input", "Input file");
  ap.parse();
  
  if (ap.files().empty())
  {
    std::cerr << "Error: input file missing" << std::endl;
    return 1;
  }
  
  if (k < 0)
  {
    std::cerr << "Error: k must be nonnegative" << std::endl;
    return 1;
  }
  
  Matrix A;
  std::ifstream inA(ap.files()[0]);
  if (!inA.good())
  {
    std::cerr << "Error: could not open '" << ap.files()[0] << "' for reading" << std::endl;
    return 1;
  }
  inA >> A;
  inA.close();
  
  Matrix::ViolationList list;
  A.identifyViolations(k, list);
  
  IntPairSet violationEntries;
  for (const Matrix::Violation& violation : list)
  {
    std::cout << "Condition " << violation._condition
              << " ; "<< "p = " << violation._p << " ; q = " << violation._q
              << " ; r = " << violation._r << " ; c = " << violation._c
              << " ; d = " << violation._d << std::endl;
    std::cout << A.getEntry(violation._p, violation._c)
              << " " << A.getEntry(violation._p, violation._d) << std::endl;
    std::cout << A.getEntry(violation._q, violation._c)
              << " " << A.getEntry(violation._q, violation._d) << std::endl;
    std::cout << A.getEntry(violation._r, violation._c)
              << " " << A.getEntry(violation._r, violation._d) << std::endl;
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
            << " / " << A.getNrTaxa() * A.getNrCharacters()
            << " = " << (double) violationEntries.size() / (A.getNrTaxa() * A.getNrCharacters()) << std::endl;
  
  return 0;
}
