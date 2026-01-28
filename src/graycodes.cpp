// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#include <bitset>
#include <limits>
#include <cassert>

#include "graycodes.h"

using namespace std;

GrayCodes::GrayCodes(int l) {
	assert(l <= numeric_limits<uint32_t>::digits);
	this->length = l;
	this->s = ~((uint32_t)0);
	this->c = 0;
	this->i = -1;
	this->changed_bit = -1;
}


bool GrayCodes::has_next() {
	return i < length;
}


uint32_t GrayCodes::get_next(int* changed_bit) {
	uint32_t result = c;
	if (changed_bit != 0) {
		*changed_bit = this->changed_bit;
	}
	i = 0;
	while (i < this->length) {
		uint32_t mask = ((uint32_t)1) << i;
		if (((c&mask) ^ (s&mask)) != 0) {
			c = c ^ mask;
			this->changed_bit = i;
			break;
		}
		s = s ^ mask;
		i += 1;
	}
	return result;
}

vector<bool> GrayCodes::to_binary(uint32_t n) {
	/**
	 * This returns the binary representation of the current state of the Gray Code ordering.
	 * index 0 contains bip for the first read.
	 */
    vector<bool> binary_vector;
	binary_vector.resize(length);
	for (uint32_t index = 0; index < length; index++) {
        // Check LSB using bitwise AND
        binary_vector[index] = ((n & 1) == 1);
        // Shift right using bitwise shift
        n >>= 1;
    }
	return binary_vector;
}