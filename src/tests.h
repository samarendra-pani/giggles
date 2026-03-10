/**
 * @file run_tests.h
 * @brief Header file for running unit tests for the C++ code.
 * @details This file calls the test files writen under src/tests/ to perform unit testing of various components of the project.
 * It includes test cases for classes and functions defined in the project to ensure correctness and reliability.
 */

#ifndef TESTS_H
#define TESTS_H

#include "tests/test_binomial.h"
#include "tests/test_column.h"
#include "tests/test_bipartitioniterator.h"
#include "tests/test_columniterator.h"
#include "tests/test_entry.h"
#include "tests/test_genotype.h"
// #include "tests/test_genotypelikelihoods.h"
#include "tests/test_graycodes.h"
#include "tests/test_variantinfo.h"
// #include "tests/test_read.h"
// #include "tests/test_readset.h"
// #include "tests/test_transitionprobabilitycomputer.h"

#include "phasing/readbipartitioning/tests/test_componentfinder.h"
#include "phasing/readbipartitioning/tests/test_phasesetcomputer.h"

#include "phasing/tests/test_phasingcolumniterator.h"
#include "phasing/tests/test_phasingcolumncostcomputer.h"


void run_all_tests() {
    std::cout << "\n================ Testing Binomial ==============================" << std::endl;
    test_binomial();
    std::cout << "\n================ Testing GrayCodes =============================" << std::endl;
    test_graycodes();
    std::cout << "\n================ Testing VariantInfo ===========================" << std::endl;
    test_variantinfo();
    std::cout << "\n================ Testing Entry =================================" << std::endl;
    test_entry();
    std::cout << "\n================ Testing Genotype ==============================" << std::endl;
    test_genotype();
    std::cout << "\n================ Testing Column ================================" << std::endl;
    test_column();
    std::cout << "\n================ Testing ColumnIterator ========================" << std::endl;
    test_columniterator();
    std::cout << "\n================ Testing BipartitionIterator ===================" << std::endl;
    test_bipartitioniterator();
    
//    test_genotypelikelihoods();
//    test_read();
//    test_readset();
//    test_transitionprobabilitycomputer();

    std::cout << "\n================ Testing ComponentFinder =======================" << std::endl;
    test_componentfinder();
//    std::cout << "\n================ Testing PhasesetComputer ======================" << std::endl;
//    test_phasesetcomputer();
    std::cout << "\n================ Testing PhasingColumnIterator =================" << std::endl;
    test_phasingcolumniterator();
    std::cout << "\n================ Testing PhasingColumnCostComputer =============" << std::endl;
    test_phasingcolumncostcomputer();
}

#endif // TESTS_H