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
int WFAWrapper::align(const std::string& text, const std::string& pattern, uint8_t type) {
    WFAligner::AlignmentStatus status;
    int score;
    int text_length;
    switch (type) {
        case 0:
            /**
             * full pattern is aligned with full text
             */
            status = aligner->alignEnd2End(pattern, text);
            score = aligner->getAlignmentScore();
            break;

        case 1:
            /**
             * full pattern is aligned to beginning part of the text
             */
            text_length = (int)text.size();
            status = aligner->alignEndsFree(pattern, 0, 0, text, 0, text_length);
            score = aligner->getAlignmentScore();
            break;
            
        case 2:
            /**
             * full pattern is aligned to end part of the text
             * 
             * Since we can start anywhere on the text but end at the end of the text,
             * the band needs to be shifted so that it can end at the end of the text.
             */
            text_length = (int)text.size();
            status = aligner->alignEndsFree(pattern, 0, 0, text, text_length, 0);
            score = aligner->getAlignmentScore();
            break;
            
        case 3:
            /**
             * full pattern is aligned to part of the text
             */
            text_length = (int)text.size();
            status = aligner->alignEndsFree(pattern, 0, 0, text, text_length, text_length);
            score = aligner->getAlignmentScore();
            break;

        default:
            break;
    }
    if ((int)status < 0) {
        throw std::runtime_error("Error: Non-success status for WFAWrapper.");
    }
    return score;
}

WFAWrapper::~WFAWrapper() {
    delete aligner;
}