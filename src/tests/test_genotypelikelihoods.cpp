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

void test_genotypelikelihood_phredscores() {
    GenotypeLikelihoods gl = GenotypeLikelihoods(4, 2);
    std::vector<uint32_t> phred_scores = gl.getPhredScores();
    for (uint32_t s: phred_scores) {
        assert(s == 0);
    }
    assert_msg(true, "GenotypeLikelihoods", "Phred scores of 0-initialized object.");
    
    std::vector<long double> likelihoods = {0.05L, 0.1L, 0.05L, 0.2L, 0.1L, 0.3L, 0.05L, 0.05L, 0.05L, 0.05L};
    gl = GenotypeLikelihoods(likelihoods, 4, 2);
    uint32_t score;

    score = gl.getPhredScore(Genotype(0, 2));
    assert(score == 7);
    score = gl.getPhredScore(Genotype(1, 2));
    assert(score == 4);
    score = gl.getPhredScore(Genotype(2, 2));
    assert(score == 7);
    score = gl.getPhredScore(Genotype(3, 2));
    assert(score == 1);
    score = gl.getPhredScore(Genotype(4, 2));
    assert(score == 4);
    score = gl.getPhredScore(Genotype(5, 2));
    assert(score == 0);
    score = gl.getPhredScore(Genotype(6, 2));
    assert(score == 7);
    score = gl.getPhredScore(Genotype(7, 2));
    assert(score == 7);
    score = gl.getPhredScore(Genotype(8, 2));
    assert(score == 7);
    score = gl.getPhredScore(Genotype(9, 2));
    assert(score == 7);

    assert_msg(true, "GenotypeLikelihoods", "Phred scores.");
}

void test_genotypelikelihoods() {

    test_genotypelikelihood_size();
    test_genotypelikelihood_phredscores();
}