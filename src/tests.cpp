#include "tests.h"

void run_all_tests() {
    std::cout << "\n================ Testing Binomial ==============================" << std::endl;
    test_binomial();
    std::cout << "\n================ Testing GrayCodes =============================" << std::endl;
    test_graycodes();
    std::cout << "\n================ Testing VariantInfo ===========================" << std::endl;
    test_variantinfo();
    std::cout << "\n================ Testing Entry =================================" << std::endl;
    test_entry();
    std::cout << "\n================ Testing Genotype ==============================" << std::endl;
    test_genotype();
    std::cout << "\n================ Testing GenotypeLikelihoods ===================" << std::endl;
    test_genotypelikelihoods();
    std::cout << "\n================ Testing Column ================================" << std::endl;
    test_column();
    std::cout << "\n================ Testing ColumnIterator ========================" << std::endl;
    test_columniterator();
    std::cout << "\n================ Testing BipartitionIterator ===================" << std::endl;
    test_bipartitioniterator();
    std::cout << "\n================ Testing HaplotypeMapper =======================" << std::endl;
    test_haplotypemapper();

//    test_read();
//    test_readset();
//    test_transitionprobabilitycomputer();

    std::cout << "\n================ Testing ComponentFinder =======================" << std::endl;
    test_componentfinder();
//    std::cout << "\n================ Testing PhasesetComputer ======================" << std::endl;
//    test_phasesetcomputer();
    std::cout << "\n================ Testing PhasingColumnIterator =================" << std::endl;
    test_phasingcolumniterator();
    std::cout << "\n================ Testing PhasingColumnCostComputer =============" << std::endl;
    test_phasingcolumncostcomputer();
}