#ifndef TEST_GRAYCODES_H
#define TEST_GRAYCODES_H

/**
 * Unit tests for GrayCodes utility functions.
 */

#include "../graycodes.h"
#include "../tests_data.h"
#include <cassert>

void test_zero_length_graycode() {
    std::cout << "[GrayCodes] Testing zero-length Gray code..." << std::endl;
    GrayCodes gray(0);
    assert_msg(gray.has_next(), "GrayCodes", "Zero-length Gray code should have one code (0).");
    assert_msg(gray.get_next() == 0, "GrayCodes", "Zero-length Gray code returns 0.");
    assert_msg(!gray.has_next(), "GrayCodes", "Zero-length Gray code has no code left..");
}

void test_l_length_graycode(int l) {
    std::cout << "[GrayCodes] Testing lengths for 1 to" << l << " Gray code..." << std::endl;
    GrayCodes gray(l);
    uint32_t expected_num_codes = 1 << l; // 2^l
    for (uint32_t i = 0; i < expected_num_codes; i++) {
        assert(gray.has_next());
        int bit_changed = -1;
        uint32_t code = gray.get_next(&bit_changed);
    }
    assert(!gray.has_next());
    assert_msg(true, "GrayCodes", "Graycode length " + std::to_string(l) + " test passed.");
}

void test_length_4_graycode() {
    std::cout << "[GrayCodes] Testing length 4 Gray code with expected sequence..." << std::endl;
    test_l_length_graycode(4);
    GrayCodes gray(4);
    std::vector<uint32_t> expected_codes = {0, 1, 3, 2, 6, 7, 5, 4, 12, 13, 15, 14, 10, 11, 9, 8};
    std::vector<int> expected_changed_bits = {-1, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0};
    uint32_t count = 0;
    while (gray.has_next()) {
        int bit_changed = -1;
        uint32_t code = gray.get_next(&bit_changed);
        assert(bit_changed == expected_changed_bits[count]);
        assert(code == expected_codes[count]);
        count += 1;
    }
    assert(count == expected_codes.size());
    assert_msg(true, "GrayCodes", "Length 4 Graycode sequence test passed.");
}

void test_graycodes() {
    test_zero_length_graycode();
    for (int l = 1; l <= 10; l++) {
        test_l_length_graycode(l);
    }
    test_length_4_graycode();
}


#endif // TEST_GRAYCODES_H