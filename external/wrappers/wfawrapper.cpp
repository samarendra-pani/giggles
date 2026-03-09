#include "wfawrapper.h"

WFAWrapper::WFAWrapper(uint32_t bandwidth) {
    aligner = new WFAlignerEdit(WFAligner::AlignmentScope::Score, WFAligner::MemoryModel::MemoryUltralow);
    aligner->setHeuristicBandedStatic(bandwidth, bandwidth);
}

/**
 * Pattern: read segment
 * Text: allele
 */
uint32_t WFAWrapper::align(const std::string& text, const std::string& pattern, uint8_t type) {
    WFAligner::AlignmentStatus status;
    uint32_t score;
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
            int text_length = (int)text.size();
            status = aligner->alignEndsFree(pattern, 0, 0, text, 0, text_length);
            score = aligner->getAlignmentScore();
            break;
            
        case 2:
            /**
             * full pattern is aligned to end part of the text
             */
            int text_length = (int)text.size();
            status = aligner->alignEndsFree(pattern, 0, 0, text, text_length, 0);
            score = aligner->getAlignmentScore();
            break;
            
        case 3:
            /**
             * full pattern is aligned to part of the text
             */
            int text_length = (int)text.size();
            status = aligner->alignEndsFree(pattern, 0, 0, text, text_length, text_length);
            score = aligner->getAlignmentScore();
            break;

        default:
            break;
    }

    return score;
}

WFAWrapper::~WFAWrapper() {
    delete aligner;
}