/**
 * @file run_tests.h
 * @brief Header file for running unit tests for the C++ code.
 * @details This file calls the test files writen under src/tests/ to perform unit testing of various components of the project.
 * It includes test cases for classes and functions defined in the project to ensure correctness and reliability.
 */

#ifndef TESTS_H
#define TESTS_H

#include "tests/test_backwardcolumniterator.h"
#include "tests/test_binomial.h"
#include "tests/test_column.h"
#include "tests/test_bipartitioniterator.h"
#include "tests/test_columniterator.h"
#include "tests/test_entry.h"
#include "tests/test_genotype.h"
#include "tests/test_genotypelikelihoods.h"
#include "tests/test_graycodes.h"
#include "tests/test_read.h"
#include "tests/test_readset.h"
#include "tests/test_transitionprobabilitycomputer.h"

#include "phasing/tests/test_phasingcolumniterator.h"


void run_all_tests() {
    test_backwardcolumniterator();
    test_binomial();
    test_column();
    test_bipartitioniterator();
    test_columniterator();
    test_entry();
//    test_genotype();
//    test_genotypelikelihoods();
    test_graycodes();
//    test_read();
//    test_readset();
    test_transitionprobabilitycomputer();

    test_phasingcolumniterator();
}

#endif // TESTS_H