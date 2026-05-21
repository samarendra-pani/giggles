#include "wfawrapper.h"
#include <iostream>
#include <stdexcept>

WFAWrapper::WFAWrapper(int32_t bandwidth): bandwidth(bandwidth) {
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
    std::string raw_cigar;
    switch (type) {
        case 0:
            /**
             * full pattern is aligned with full text
             */
            if (abs((int)text.size() - (int)pattern.size()) > bandwidth) {
                throw std::runtime_error("Error: Difference in text and pattern length is greater than bandwidth.");
            }
            aligner->setHeuristicBandedStatic(-bandwidth, bandwidth);
            status = aligner->alignEnd2End(pattern, text);
            score = aligner->getAlignmentScore();
            aligner->setHeuristicNone();
            break;

        case 1:
            /**
             * full pattern is aligned to beginning part of the text
             */
            text_length = (int)text.size();
            aligner->setHeuristicBandedStatic(-bandwidth, bandwidth);
            status = aligner->alignEndsFree(pattern, 0, 0, text, 0, text_length);
            score = aligner->getAlignmentScore();
            aligner->setHeuristicNone();
            break;
            
        case 2:
            /**
             * full pattern is aligned to end part of the text
             * 
             * Since we can start anywhere on the text but end at the end of the text,
             * the band needs to be shifted so that it can end at the end of the text.
             */
            text_length = (int)text.size();
            aligner->setHeuristicBandedStatic(text_length-(int)pattern.size()-bandwidth, text_length-(int)pattern.size()+bandwidth);
            status = aligner->alignEndsFree(pattern, 0, 0, text, text_length, 0);
            score = aligner->getAlignmentScore();
            aligner->setHeuristicNone();
            break;
            
        case 3:
            /**
             * full pattern is aligned to part of the text
             */
            text_length = (int)text.size();
            aligner->setHeuristicNone();
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