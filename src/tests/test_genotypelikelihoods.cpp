#include "test_genotypelikelihoods.h"

void test_genotypelikelihood_size() {
    GenotypeLikelihoods gl = GenotypeLikelihoods(2, 2);
    assert(gl.size() ==  3);
    gl = GenotypeLikelihoods(3, 2);
    assert(gl.size() ==  6);
    gl = GenotypeLikelihoods(4, 2);
    assert(gl.size() ==  10);
    gl = GenotypeLikelihoods(5, 2);
    assert(gl.size() ==  15);
    gl = GenotypeLikelihoods(2, 1);
    assert(gl.size() ==  2);
    gl = GenotypeLikelihoods(3, 1);
    assert(gl.size() ==  3);
    gl = GenotypeLikelihoods(4, 1);
    assert(gl.size() ==  4);
    gl = GenotypeLikelihoods(5, 1);
    assert(gl.size() ==  5);

    assert_msg(true, "GenotypeLikelihoods", "Likelihood vector size.");
}

void test_genotypelikelihoods() {

    test_genotypelikelihood_size();
}