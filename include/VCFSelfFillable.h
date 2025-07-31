#ifndef __VCFSELFFILLABLE
#define __VCFSELFFILLABLE

#include "SelfRecordSet.h"

class VCFSelfHetero;
class Option;


//////////////////// VCFSelfFillable ////////////////////

class VCFSelfFillable : public VCFBase, public VCFSmallBase {
	using PosWithChr = std::tuple<int,ll,std::string>;
	
public:
	using Position = std::tuple<int, ll, std::string>;
	using Item = std::pair<std::vector<VCFSelfHetero *>,
							std::vector<VCFImpSelfRecord *>>;
	
	struct ConfigThreadPhase {
		const std::size_t	first;
		const std::size_t	num_threads;
		const std::vector<SelfRecordSet *>&	record_sets;
		VCFSelfFillable	*vcf;
		
		ConfigThreadPhase(std::size_t f, std::size_t n,
							const std::vector<SelfRecordSet *>& r,
							VCFSelfFillable *v) :
						first(f), num_threads(n), record_sets(r), vcf(v) { }
	};
	
private:
	std::vector<VCFSelfFillableRecord *>	records;
	
public:
	VCFSelfFillable(const std::vector<STRVEC>& h, const STRVEC& s,
							const std::vector<VCFSelfFillableRecord *>& rs);
	virtual ~VCFSelfFillable();
	
	///// virtual methods for VCFSmallBase /////
	std::size_t size() const override { return records.size(); }
	VCFRecord *get_record(std::size_t i) const override {
		return records[i];
	}
	
	///// virtual methods for VCFFamilyBase /////
	const std::vector<STRVEC>& get_header() const override {
		return VCFBase::get_header();
	}
	const STRVEC& get_samples() const override {
		return VCFBase::get_samples();
	}
	
	///// non-virtual methods /////
	VCFSelfFillableRecord *get_fillable_record(std::size_t i) const {
		return records[i];
	}
	const std::vector<VCFSelfFillableRecord *>& get_records() const {
		return records;
	}
	
	void modify();
	void set_records(const std::vector<VCFSelfFillableRecord *>& rs) {
		records = rs;
	}
	void clear_records() { records.clear(); }
	
protected:
	template<typename Iter>
	VCFSelfFillableRecord *find_neighbor_same_type_record(
					std::size_t i, std::size_t c, Iter first, Iter last) const;
	VCFSelfFillableRecord *find_prev_same_type_record(std::size_t i,
													std::size_t c) const;
	VCFSelfFillableRecord *find_next_same_type_record(std::size_t i,
													std::size_t c) const;
	
public:
	static VCFSmall *merge(const std::vector<VCFSelfFillable *>& vcfs,
											const STRVEC& orig_samples);
	static VCFSelfFillable *fill(const std::vector<VCFSelfHetero *>& vcfs,
							 const std::vector<VCFImpSelfRecord *>& records);
	static std::vector<VCFSelfFillableRecord *> merge_records(
							const std::vector<VCFSelfHetero *>& vcfs,
							const std::vector<VCFImpSelfRecord *>& records);
	
protected:
	static std::vector<std::vector<VCFSelfFillableRecord *>>
				collect_records(const std::vector<VCFSelfFillable *>& vcfs);
};
#endif
