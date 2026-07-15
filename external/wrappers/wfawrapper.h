#ifndef WFAWRAPPER_H
#define WFAWRAPPER_H

#include "../external/wfa2/bindings/cpp/WFAligner.hpp"

using namespace wfa;

/**
 * Wrapper class for WFA2-lib (https://github.com/smarco/WFA2-lib/).
 */
class WFAWrapper {
    
    public:
        
        WFAWrapper();

        virtual ~WFAWrapper();
       
        /**
         * Align allele and query.
         * Variable type specifies what type of alignment to perform:
         * 0 - end-to-end alignment
         * 1 - start-fixed alignment
         * 2 - end-fixed alignment
         * 3 - free alignment
         * 
         * @returns alignment score
         */
        int align(const char* allele, int allele_len, const char* query, int query_len, uint8_t type);

    private:
        /**
         * aligners defined based on their type
         */
        WFAlignerEdit* aligner;

        /**
         * Internal switch case logic
         * Returns negative status if status is negative. Otherwise returns score.
         */
        int attempt_alignment(const char* allele, int allele_len, const char* query, int query_len, uint8_t type);

};

#endif // WFAWRAPPER_H