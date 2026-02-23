#include "wfawrapper.h"

WFAWrapper::WFAWrapper(uint32_t bandwidth) {
    aligner0 = new WFAlignerEdit(WFAligner::AlignmentScope::Score, WFAligner::MemoryModel::MemoryUltralow);
    aligner0->setHeuristicBandedStatic(bandwidth, bandwidth);
}

uint32_t WFAWrapper::align(const std::string& text, const std::string& pattern, uint8_t type) {
    WFAligner::AlignmentStatus status;
    uint32_t score;
    switch (type)
    {
    case 0:
        status = aligner0->alignEnd2End(pattern, text);
        score = aligner0->getAlignmentScore();
        break;

    case 1:
        break;
        
    case 2:
        break;
        
    case 3:
        //status = aligner0->alignEndsFree(pattern, );
        //score = aligner0->getAlignmentScore();
        break;

    default:
        break;
    }

    return score;
}

WFAWrapper::~WFAWrapper() {
    delete aligner0;
    delete aligner1;
    delete aligner2;
    delete aligner3;
}