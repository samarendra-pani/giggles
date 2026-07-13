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
int WFAWrapper::align(const char* text, int text_len, const char* pattern, int pattern_len, uint8_t type) {
    
    int score = attempt_alignment(text, text_len, pattern, pattern_len, type);
    if ((int)score < 0) {
        // If it's a memory/step capacity issue, clear the slate and try ONE more time
        if ((int)score == -100) { 
            std::cerr << "[Core::WFAWrapper] Warning: Max steps reached. Resetting allocator and retrying..." << std::endl;
            // Recreate the aligner inside the wrapper
            delete aligner;
            aligner = new WFAlignerEdit(WFAligner::AlignmentScope::Score, WFAligner::MemoryModel::MemoryHigh);
            
            // Re-run the exact same alignment logic
            score = attempt_alignment(text, text_len, pattern, pattern_len, type); 
            if ((int)score >= 0) {
                return score;
            }
        }
        
        // If it still fails, or it's a different error, throw
        throw std::runtime_error("Error: Non-success status for WFAWrapper after retry.");
    }
    return score;
}

int WFAWrapper::attempt_alignment(const char* text, int text_len, const char* pattern, int pattern_len, uint8_t type) {
    WFAligner::AlignmentStatus status;
    int score;
    switch (type) {
        case 0:
            /**
             * full pattern is aligned with full text
             */
            status = aligner->alignEnd2End(pattern, pattern_len, text, text_len);
            score = aligner->getAlignmentScore();
            break;

        case 1:
            /**
             * full pattern is aligned to beginning part of the text
             */
            
            status = aligner->alignEndsFree(pattern, pattern_len, 0, 0, text, text_len, 0, text_len);
            score = aligner->getAlignmentScore();
            break;
            
        case 2:
            /**
             * full pattern is aligned to end part of the text
             * 
             * Since we can start anywhere on the text but end at the end of the text,
             * the band needs to be shifted so that it can end at the end of the text.
             */
            status = aligner->alignEndsFree(pattern, pattern_len, 0, 0, text, text_len, text_len, 0);
            score = aligner->getAlignmentScore();
            break;
            
        case 3:
            /**
             * full pattern is aligned to part of the text
             */
            status = aligner->alignEndsFree(pattern, pattern_len, 0, 0, text, text_len, text_len, text_len);
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