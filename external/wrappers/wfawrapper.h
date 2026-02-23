#ifndef WFAWRAPPER_H
#define WFAWRAPPER_H

#include "../external/wfa2/bindings/cpp/WFAligner.hpp"

using namespace wfa;

/**
 * Wrapper class for WFA2-lib (https://github.com/smarco/WFA2-lib/).
 * 
 * @note pywfa (https://github.com/kcleal/pywfa/) exists but does not support edit distance metric
 *   and one-end-free alignment.
 */
class WFAWrapper {
    
    public:
        
        WFAWrapper(uint32_t bandwidth);

        virtual ~WFAWrapper();
       
        /**
         * Align text and pattern.
         * Variable type specifies what type of alignment to perform:
         * 0 - end-to-end alignment
         * 1 - start-fixed alignment
         * 2 - end-fixed alignment
         * 3 - free alignment
         * 
         * @returns alignment score
         */
        uint32_t align(const std::string& text, const std::string& pattern, uint8_t type);

    private:
        /**
         * aligners defined based on their type
         */
        WFAlignerEdit* aligner0;
        WFAlignerEdit* aligner1;
        WFAlignerEdit* aligner2;
        WFAlignerEdit* aligner3;

};

#endif // WFAWRAPPER_H