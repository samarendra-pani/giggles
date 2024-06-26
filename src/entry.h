// Code modified from WhatsHap (https://github.com/whatshap/whatshap)

#ifndef ENTRY_H
#define ENTRY_H

#include <iostream>
#include<vector>

class Entry {
public:
	
	Entry(unsigned int r, int m, std::vector<double> e, int q, int reg_const, double base_const);
	
	// For blank entries.
	Entry(unsigned int r, int m);

	unsigned int get_read_id() const;
	int get_allele_type() const;
	std::vector<long double> get_emission_score() const;
	int get_quality() const;

	void set_read_id(unsigned int r);
	void set_allele_type(int m);
	void set_emission_score(std::vector<double> e, int reg_const, double base_const);
	void set_quality(int q);

	friend std::ostream& operator<<(std::ostream& out, const Entry& e);

private:
	unsigned int read_id;
	int allele;
	std::vector<long double> emission_score;
	int quality;
};

#endif
