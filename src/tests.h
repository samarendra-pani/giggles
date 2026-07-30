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
#include "tests/test_genotypelikelihoods.h"
#include "tests/test_graycodes.h"
#include "tests/test_haplotypemapper.h"
#include "tests/test_variantinfo.h"
// #include "tests/test_read.h"
#include "tests/test_readset.h"
// #include "tests/test_transitionprobabilitycomputer.h"

#include "phasing/readbipartitioning/tests/test_componentfinder.h"
#include "phasing/readbipartitioning/tests/test_phasesetcomputer.h"
#include "phasing/readbipartitioning/tests/test_haplotagcomputer.h"
#include "phasing/readbipartitioning/tests/test_set_cluster_ids.h"

#include "phasing/tests/test_phasingcolumniterator.h"
#include "phasing/tests/test_phasingcolumncostcomputer.h"
#include "phasing/tests/test_phasingdptable.h"

#include "haplotypesampler/tests/test_haplotypesampler.h"

void run_all_tests();

#endif // TESTS_H