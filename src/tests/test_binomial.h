#ifndef TEST_BINOMIAL_H
#define TEST_BINOMIAL_H

/**
 * Unit tests for binomial functions.
 */

#include "../binomial.h"
#include "../tests_data.h"
#include <cassert>

void test_binomial() {
    // test binomial coefficient
    assert(binomial_coefficient(5, 2) == 10);
    assert(binomial_coefficient(10, 3) == 120);
    assert(binomial_coefficient(0, 0) == 1);
    assert(binomial_coefficient(7, 0) == 1);
    assert(binomial_coefficient(7, 7) == 1);
    assert(binomial_coefficient(6, 4) == 15);

    assert_msg(true, "Binomial", "All binomial coefficient tests passed.");
}

#endif // TEST_BINOMIAL_H