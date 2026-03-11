#ifndef TESTS_DATA_H
#define TESTS_DATA_H

#include <cassert>

#include "variantinfo.h"
#include "readset.h"

// Helper to print success messages
void assert_msg(bool condition, const std::string& prefix, const std::string& message);

std::vector<variant_information_t>  mock_variant_info_table_1();

ReadSet* mock_readset_1();

std::vector<variant_information_t> mock_variant_info_table_2();

ReadSet* mock_readset_2();

ReadSet* mock_superreads();

#endif // TESTS_DATA_H