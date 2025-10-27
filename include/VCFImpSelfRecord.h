#ifndef __VCFIMPSELFRECORD
#define __VCFIMPSELFRECORD

#include <vector>

#include "GenoRecord.h"
#include "ClassifyRecord.h"


//////////////////// SelfFillType ////////////////////

enum class SelfFillType {
	P01, FILLED, IMPUTABLE, UNABLE
};


//////////////////// VCFImpSelfRecord ////////////////////

class VCFImpSelfRecord : public GenoRecord {
public:
	typedef std::pair<VCFImpSelfRecord *, int>	RecordWithPos;
	
protected:
	const int	index;
	WrongType	wrong_type;
	const ParentComb	comb;
	
public:
	VCFImpSelfRecord(ll pos, const std::vector<int>& geno,
							int i, WrongType type, ParentComb c) :
									GenoRecord(pos, geno), index(i),
									wrong_type(type), comb(c) { }
	virtual ~VCFImpSelfRecord() { }
	
	int get_index() const { return index; }
	ParentComb get_comb() const { return comb; }
	
	bool is_right() const {
		return this->wrong_type == WrongType::RIGHT;
	}
	bool is_fixed() const {
		return (this->comb == ParentComb::P00x00 ||
					this->comb == ParentComb::P11x11) ||
				((this->comb == ParentComb::P00x01 ||
					this->comb == ParentComb::P01x11) && is_right());
	}
	
	virtual bool is_imputable() const = 0;
	virtual SelfFillType get_fill_type() const = 0;
	
	void enable_modification() {
		this->wrong_type = WrongType::MODIFIABLE;
	}
	
public:
	static int get_records_type(const std::vector<RecordWithPos>& records);
	static std::size_t find_fixed_index(
		const std::vector<std::pair<int, std::vector<RecordWithPos>>>& items);
	static std::pair<int, std::vector<RecordWithPos>>
	which_is_fixed(const std::vector<RecordWithPos>& v);
};
#endif
