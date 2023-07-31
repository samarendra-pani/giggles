#include <cassert>
#include <bits/stdc++.h>
#include <cmath>
using namespace std;

#include "entry.h"

Entry::Entry(unsigned int r, int m, std::vector<double> e, int q) {
	read_id = r;
	allele = m;
	set_emission_score(e);
	quality = q;
}

Entry::Entry(unsigned int r, int m) {
	read_id = r;
	allele = m;
}

unsigned int Entry::get_read_id() const {
	return read_id;
}


int Entry::get_allele_type() const {
	return allele;
}


std::vector<long double> Entry::get_emission_score() const {
	return emission_score;
}

int Entry::get_quality() const {
	return quality;
}

void Entry::set_read_id(unsigned int r) {
	read_id = r;
}


void Entry::set_allele_type(int m) {
	allele = m;
}


void Entry::set_emission_score(std::vector<double> e) {
	emission_score.resize(e.size());
	int i = 0;
	double normalization = 0.0L;
	/* for (auto it = e.begin(); it != e.end(); it++, i++) {
		emission_score[i] = 1e-10L;
		long double score = 0.0L;
		for (auto it2 = e.begin(); it2 != e.end(); it2++) {
			score += exp(*it2 - *it);
		}
		emission_score[i] += 1/score;
		normalization += emission_score[i];
	} */
	for (auto it = e.begin(); it != e.end(); it++, i++) {
		emission_score[i] = 1e-10L;
		emission_score[i] += 1/(-(*it)+1);
		normalization += emission_score[i];
	}
	/* cout << "Normalization: " << normalization << endl;;
	for (auto it: emission_score){
		cout << it << ' ';
	}
	cout << endl; */
	transform((emission_score).begin(), (emission_score).end(), (emission_score).begin(), std::bind2nd(std::divides<long double>(), normalization));
	/* for (auto it: e){
		cout << it << ' ';
	}
	cout << endl;
	for (auto it: emission_score){
		cout << it << ' ';
	}
	cout << endl;
	cout << endl;
	auto max = *max_element(emission_score.begin(), emission_score.end());
	cout << max << endl;
	if (max != emission_score[0]) exit(0); */
}
		

void Entry::set_quality(int q) {
	quality = q;
}

std::ostream& operator<<(std::ostream& out, const Entry& e) {
	out << "Entry(" << e.read_id ;
	out << ","<< e.allele << ",(";
	for (auto i : e.emission_score) {
		out << i << ",";
	}
	out << ")," << e.quality << ")" << std::endl;
	return out;
}
