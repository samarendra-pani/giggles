#ifndef WFAWRAPPER_H
#define WFAWRAPPER_H

#include "../external/wfa2/bindings/cpp/WFAligner.hpp"

using namespace wfa;

/**
 * Wrapper class for WFA2-lib (https://github.com/smarco/WFA2-lib/).
 */
class WFAWrapper {
    
    public:
        
        WFAWrapper(int32_t bandwidth);

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
        int align(const std::string& text, const std::string& pattern, uint8_t type);

    private:
        /**
         * aligners defined based on their type
         */
        WFAlignerEdit* aligner;

        /** the bandwidth */
        int32_t bandwidth;

};

#endif // WFAWRAPPER_H