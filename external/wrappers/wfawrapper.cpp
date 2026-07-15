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
int WFAWrapper::align(const char* allele, int allele_len, const char* query, int query_len, uint8_t type) {
    
    int score = attempt_alignment(allele, allele_len, query, query_len, type);
    if ((int)score < 0) {
        // If it's a memory/step capacity issue, clear the slate and try ONE more time
        if ((int)score == -100) { 
            std::cerr << "[Core::WFAWrapper] Warning: Max steps reached. Resetting allocator and retrying..." << std::endl;
            // Recreate the aligner inside the wrapper
            delete aligner;
            aligner = new WFAlignerEdit(WFAligner::AlignmentScope::Score, WFAligner::MemoryModel::MemoryHigh);
            
            // Re-run the exact same alignment logic
            score = attempt_alignment(allele, allele_len, query, query_len, type); 
            if ((int)score >= 0) {
                return score;
            }
        }
        
        // If it still fails, or it's a different error, throw
        std::string error = "Error: Non-success status " + std::to_string((int)score) + " for WFAWrapper after retry.";
        throw std::runtime_error(error);
    }
    return score;
}

int WFAWrapper::attempt_alignment(const char* allele, int allele_len, const char* query, int query_len, uint8_t type) {
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
            status = aligner->alignEndsFree(query, query_len, 0, 0, allele, allele_len, 0, allele_len);
            score = aligner->getAlignmentScore();
            //aligner->setHeuristicNone();
            break;

        case 3:
            /**
             * full query is aligned to part of the allele
             */
            status = aligner->alignEndsFree(query, query_len, 0, 0, allele, allele_len, allele_len, allele_len);
            score = aligner->getAlignmentScore();
            break;

        default:
            break;
    }
    if ((int)status < 0) {
        return (int)status;
    }
    return score;
}

WFAWrapper::~WFAWrapper() {
    delete aligner;
}