#include "wfawrapper.h"
#include <iostream>
#include <stdexcept>

WFAWrapper::WFAWrapper() {
    aligner = new WFAlignerEdit(WFAligner::AlignmentScope::Score, WFAligner::MemoryModel::MemoryHigh);
}

/**
 * Pattern: read segment
 * Text: allele
 */
int WFAWrapper::align(const char* query, int query_len, const char* allele, int allele_len, uint8_t type) {
    
    int score = attempt_alignment(query, query_len, allele, allele_len, type);
    if ((int)score < 0) {
        std::string error = "Error: Non-success status " + std::to_string((int)score) + " for WFAWrapper after retry.";
        throw std::runtime_error(error);
    }
    return score;
}

int WFAWrapper::attempt_alignment(const char* query, int query_len, const char* allele, int allele_len, uint8_t type) {
    WFAligner::AlignmentStatus status;
    int score;
    int bandwidth;
    switch (type) {
        case 0:
            /**
             * full query is aligned with full allele
             */
            bandwidth = (int)(std::max(allele_len, query_len)/6.0);
            aligner->setHeuristicBandedStatic(-bandwidth, bandwidth);
            aligner->setMaxAlignmentSteps(bandwidth);
            status = aligner->alignEnd2End(query, query_len, allele, allele_len);
            score = aligner->getAlignmentScore();
            //aligner->setHeuristicNone();
            break;
    
        case 1:
            /**
             * full query is aligned to end part of the allele
             * 
             * Since we can start anywhere on the allele but end at the end of the allele,
             * the band needs to be shifted so that it can end at the end of the allele.
             */
            bandwidth = (int)(query_len/6.0);
            aligner->setHeuristicBandedStatic(allele_len-query_len-bandwidth, allele_len-query_len+bandwidth);
            aligner->setMaxAlignmentSteps(bandwidth);
            status = aligner->alignEndsFree(query, query_len, 0, 0, allele, allele_len, allele_len, 0);
            score = aligner->getAlignmentScore();
            //aligner->setHeuristicNone();
            break;
        
        case 2:
            /**
             * full query is aligned to beginning part of the allele
             */
            bandwidth = (int)(query_len/6.0);
            aligner->setHeuristicBandedStatic(-bandwidth, bandwidth);
            aligner->setMaxAlignmentSteps(bandwidth);
            status = aligner->alignEndsFree(query, query_len, 0, 0, allele, allele_len, 0, allele_len);
            score = aligner->getAlignmentScore();
            //aligner->setHeuristicNone();
            break;

        case 3:
            /**
             * full query is aligned to part of the allele
             */
            bandwidth = (int)(query_len/6.0);
            aligner->setHeuristicNone();
            aligner->setMaxAlignmentSteps(bandwidth);
            status = aligner->alignEndsFree(query, query_len, 0, 0, allele, allele_len, allele_len, allele_len);
            score = aligner->getAlignmentScore();
            break;

        default:
            break;
    }
    if ((int)status == -100) {
        return bandwidth;
    }
    if ((int)status < 0) {
        return (int)status;
    }
    return score;
}

WFAWrapper::~WFAWrapper() {
    delete aligner;
}