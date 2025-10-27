#include <stdexcept>
#include "../include/Genotype.h"
#include "../include/common.h"

using namespace std;


//////////////////// Genotype ////////////////////

const int Genotype::UN_00;
const int Genotype::UN_01;
const int Genotype::UN_11;
const int Genotype::NA;
const int Genotype::PH_00;
const int Genotype::PH_01;
const int Genotype::PH_10;
const int Genotype::PH_11;

Genotype::Genotype(const string& s) : gt1(s.c_str()[0]), gt2(s.c_str()[2]),
												phasing(s.c_str()[1] == '|') { }

pair<char,char> Genotype::gts() const {
	return pair<char,char>(gt1, gt2);
}

bool Genotype::includes(char gt) const {
	return gt == gt1 || gt == gt2;
}

string Genotype::int_to_gt(int n) {
	switch(n) {
		case 0:  return "0/0";
		case 1:  return "0/1";
		case 2:  return "1/1";
		default: return "./.";
	}
}

bool Genotype::conflicts(const Genotype& mat_gt, const Genotype& pat_gt,
												bool considers_phasing) {
	if(considers_phasing && this->phasing) {
		return mat_gt.includes(gt1) && pat_gt.includes(gt2);
	}
	else {
		if(gt1 == gt2)
			return mat_gt.includes(gt1) && pat_gt.includes(gt2);
		else
			return (mat_gt.includes(gt1) && pat_gt.includes(gt2)) ||
					(pat_gt.includes(gt1) && mat_gt.includes(gt2));
	}
}

bool Genotype::is_valid(int gt, int mat_gt, int pat_gt) {
	const vector<int>	mat_gts = possible_gts(mat_gt);
	const vector<int>	pat_gts = possible_gts(pat_gt);
	for(auto p = mat_gts.begin(); p != mat_gts.end(); ++p) {
		for(auto q = pat_gts.begin(); q != pat_gts.end(); ++q) {
			if(*p + *q == 0 && Genotype::is_00(gt))
				return true;
			else if(*p + *q == 1 && Genotype::is_01(gt))
				return true;
			else if(*p + *q == 2 && Genotype::is_11(gt))
				return true;
		}
	}
	return false;
}

vector<int> Genotype::possible_gts(int gt) {
	if(Genotype::is_00(gt))
		return {0};
	else if(Genotype::is_11(gt))
		return {1};
	else
		return {0, 1};
}

size_t Genotype::find_key_position(const string& info, const string& key) {
	const auto	keys = Common::split(info, ':');
	for(size_t i = 0; i < keys.size(); ++i) {
		if(keys[i] == key) {
			return i;
		}
	}
	return string::npos;
}
