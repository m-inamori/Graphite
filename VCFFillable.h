#ifndef __VCFFILLABLE
#define __VCFFILLABLE

#include "VCFImpFamily.h"
#include "ClassifyRecord.h"

class VCFHeteroHomo;
class Option;


//////////////////// Genotype ////////////////////

class Genotype {
	const char	gt1;
	const char	gt2;
	const bool	phasing;
	
public:
	Genotype(const std::string& s);
	
	std::pair<char,char> gts() const;
	bool includes(char gt) const;
	bool conflicts(const Genotype& mat_gt, const Genotype& pat_gt,
													bool considers_phasing);
	
	static bool is_valid(const std::string& gt, int mat_gt, int pat_gt);
	static std::string possible_gts(int gt);
	static int sum_gt(const std::string& gt);
	static bool is_all_NA(const std::vector<std::string>& GTs);
};


//////////////////// VCFFillableRecord ////////////////////

class VCFFillableRecord : public VCFFamilyRecord {
protected:
	int			index;
	FillType	type;
	ParentComb	comb;
	
public:
	VCFFillableRecord(const STRVEC& v, const STRVEC& s,
									int i, FillType t, ParentComb c) :
						VCFFamilyRecord(v, s), index(i), type(t), comb(c) { }
	
	int get_index() const { return index; }
	FillType get_type() const { return type; }
	bool is_unable() const { return type == FillType::UNABLE; }
	bool is_fillable_type() const {
		return type == FillType::UNABLE || type == FillType::IMPUTABLE;
	}
	bool is_filled_type() const {
		return type == FillType::FILLED || type == FillType::MAT
										|| type == FillType::PAT;
	}
	bool is_mat_type() const { return type == FillType::MAT; }
	bool is_pat_type() const { return type == FillType::PAT; }
	
	ParentComb get_comb() const { return comb; }
	bool is_00x11() const { return comb == ParentComb::P00x11; }
	void set_00x11() { comb = ParentComb::P00x11; }
	
	std::vector<std::pair<int, int>>	possible_phasings() const;
	
	STRVEC prog_gts() const;
	VCFFillableRecord *copy() const;
	std::tuple<int,int,int> count_gts() const;
	std::string gt_from_parent(int mat_from, int pat_from) const;
	std::string gt_from_mat(int mat_from, int c) const;
	std::string gt_from_pat(int pat_from, int c) const;
	int mat_from(int c) const;
	int pat_from(int c) const;
	
	void modify_gts();
	int from_which_chrom(std::size_t i, bool is_mat) const;
	void modify_gts(const STRVEC& new_prog_gts);
	void fill_PGT();
	void set(const STRVEC& new_v, FillType new_type);
	void modify_parents_type();
	
private:
	std::vector<std::vector<double>> make_probability_table() const;
	
	void phase();
	int find_geno_type(const std::string& type) const;
	
	std::string inverse_gt(const std::string& gt, bool inv) const;
	std::string inverse_prog_gt(const std::string& gt,
									bool inv_mat, bool inv_pat) const;
	STRVEC inverse_prog_gts(const STRVEC& prog_gts,
									bool inv_mat, bool inv_pat) const;
	void inverse_parents_gts(bool inv_mat, bool inv_pat);
	bool is_same_gts(const std::string& gt1, const std::string& gt2) const;
	bool is_near_prog_gts(const STRVEC& gts) const;
	int hash(int d) const;
	std::string decide_by_majority(const std::vector<std::string>& GTs) const;
	void swap_parents(int i, const std::string& GT);
	
public:
	static VCFRecord *integrate_records(
						const std::vector<VCFFillableRecord *>& records);
	static std::string decide_duplicated_Genotype(
							const std::vector<VCFFillableRecord *>& records,
							const std::vector<std::pair<int, int>>& positions);
	static VCFFillableRecord *convert(const VCFImpFamilyRecord *record);
	static VCFRecord *merge(const std::vector<VCFFillableRecord *>& records,
														const STRVEC& samples);
	static void integrate_each_sample(
							const std::vector<VCFFillableRecord *>& records,
							const std::vector<std::pair<int, int>>& positions);
	static bool is_all_same_GT(const std::vector<VCFFillableRecord *>& records,
						const std::vector<std::pair<int, int>>& pos_samples);
	static VCFRecord *integrate(const std::vector<VCFFillableRecord *>& records,
			const std::vector<std::string>& samples,
			const std::vector<std::vector<std::pair<int, int>>>& pos_samples);
	static int from_which_chrom(const VCFFillableRecord *record,
										std::size_t i, bool is_mat);
	static int from_which_chrom_mat(const VCFFillableRecord *record,
														std::size_t i) {
		return from_which_chrom(record, i, true);
	}
	static int from_which_chrom_pat(const VCFFillableRecord *record,
														std::size_t i) {
		return from_which_chrom(record, i, false);
	}
};


//////////////////// VCFFillable ////////////////////

class VCFFillable : public VCFSmall {
	using Pair = std::pair<int,int>;
	using PosWithChr = std::tuple<int,ll,std::string>;
	
	class RecordSet {
	public:
		VCFFillableRecord	*record;
		VCFFillableRecord	*prev_mat_record;
		VCFFillableRecord	*next_mat_record;
		VCFFillableRecord	*prev_pat_record;
		VCFFillableRecord	*next_pat_record;
		
	public:
		RecordSet(VCFFillableRecord *r,
				  VCFFillableRecord *pm, VCFFillableRecord *nm,
				  VCFFillableRecord *pp, VCFFillableRecord *np) : record(r),
						  			prev_mat_record(pm), next_mat_record(nm),
						  			prev_pat_record(pp), next_pat_record(np) { }
		
		std::string gt_each(int i, const VCFFillableRecord *record) const;
		std::string gt(int i) const { return gt_each(i, record); }
		std::string prev_mat_gt(int i) const;
		std::string next_mat_gt(int i) const;
		std::string prev_pat_gt(int i) const;
		std::string next_pat_gt(int i) const;
		
		int prev_mat_from(std::size_t i) const;
		int next_mat_from(std::size_t i) const;
		int prev_pat_from(std::size_t i) const;
		int next_pat_from(std::size_t i) const;
		
		int near_mat_from(size_t i) const;
		int near_pat_from(size_t i) const;
		
		bool is_mat_prev_near() const;
		bool is_pat_prev_near() const;
		
		Pair select_nearest_froms(const std::vector<Pair>& pairs,
													std::size_t i) const;
		Pair select_pair(const std::vector<Pair>& pairs,
								std::size_t i, bool selected=false) const;
		std::vector<double> probs_from_which_chrom(int prev_chrom,
													int next_chrom) const;
		std::vector<double> probs_from_which_chrom(std::size_t i,
													bool is_mat) const;
		double compute_phasing_likelihood_each(std::size_t i,
									int mat_phasing, int pat_phasing) const;
		double likelihood_each(const std::string& gt,
								const std::vector<double>& probs_mat,
								const std::vector<double>& probs_pat,
								int mat_phasing, int pat_phasing) const;
		double compute_phasing_likelihood(int mat_phasing,
											int pat_phasing) const;
		void determine_phasing();
		int select_from(int from1, int from2,
									const VCFRecord *record1,
									const VCFRecord *record2) const;
		std::string modify_gt(size_t i);
		void impute_core();
		
		int select_mat(const std::vector<Pair>& pairs) const;
		void impute_NA_mat_each(std::size_t i) const;
		int select_pat(const std::vector<Pair>& pairs) const;
		void impute_NA_pat_each(std::size_t i) const;
	};
	
	struct ConfigReplaceThread {
		const std::vector<VCFFillable *>&	vcfs;
		const std::vector<const VCFRecord *>&	orig_records;
		const std::size_t	first;
		const int	num_thread;
		
		ConfigReplaceThread(const std::vector<VCFFillable *>& v,
						const std::vector<const VCFRecord *>& o, int f, int n) :
					vcfs(v), orig_records(o), first(f), num_thread(n) { }
	};
	
	struct ConfigThread {
		const std::vector<VCFFillable *>&	vcfs;
		const std::size_t	first;
		const int	num_thread;
		
		ConfigThread(const std::vector<VCFFillable *>& v, int f, int n) :
								vcfs(v), first(f), num_thread(n) { }
	};
	
	using Position = std::tuple<int, ll, std::string>;
	using Group = std::pair<FillType, std::vector<VCFFillableRecord *>>;
	
	std::vector<VCFFillableRecord *>	fillable_records;
	
public:
	VCFFillable(const std::vector<STRVEC>& h, const STRVEC& s,
									std::vector<VCFFillableRecord *> rs);
	~VCFFillable() { }
	
	void modify();
	VCFFillable *create_from_header() const;
	void set_records(const std::vector<VCFFillableRecord *>& rs);
	void set_records_base(const std::vector<VCFFillableRecord *>& rs);
	
	VCFFillableRecord *get_record(std::size_t i) const {
		return fillable_records[i];
	}
	const std::vector<VCFFillableRecord *>& get_records() const {
		return fillable_records;
	}

private:
	std::vector<Group> group_records() const;
	VCFFillableRecord *find_prev_record(FillType type, int i,
										const std::vector<Group>& groups) const;
	VCFFillableRecord *find_next_record(FillType type, int i,
										const std::vector<Group>& groups) const;
	void phase(int i, const std::vector<Group>& groups);
	
	template<typename Iter>
	VCFFillableRecord *find_neighbor_same_type_record(
					std::size_t i, std::size_t c, Iter first, Iter last) const;
	VCFFillableRecord *find_prev_same_type_record(std::size_t i,
													std::size_t c) const;
	VCFFillableRecord *find_next_same_type_record(std::size_t i,
													std::size_t c) const;
	
	const RecordSet *create_recordset(std::size_t i,
										std::size_t c, bool is_mat) const;
	void impute_NA_mat_each(std::size_t i, std::size_t c);
	void impute_NA_mat(std::size_t i);
	
	void impute_NA_pat_each(std::size_t i, std::size_t c);
	void impute_NA_pat(std::size_t i);
	
	std::pair<int, int> find_prev_mat_from(int i, int c) const;
	std::pair<int, int> find_next_mat_from(int i, int c) const;
	std::pair<int, int> find_prev_pat_from(int i, int c) const;
	std::pair<int, int> find_next_pat_from(int i, int c) const;
	int select_from(const std::pair<int, int>& f1,
					const std::pair<int, int>& f2, int i) const;
	int find_mat_from(int i, int c) const;
	int find_pat_from(int i, int c) const;
	void impute_others(int i);
	
public:
	static VCFFillable *fill(const std::vector<VCFHeteroHomo *>& vcfs,
				const std::vector<VCFImpFamilyRecord *>& records, bool all_out);
	static std::vector<VCFFillableRecord *> merge_records(
							const std::vector<VCFHeteroHomo *>& vcfs,
							const std::vector<VCFImpFamilyRecord *>& records,
							bool all_out);
	static std::vector<std::vector<VCFFillableRecord *>>
		collect_records(const std::vector<VCFFillable *>& vcfs, bool all_out);
	static std::pair<STRVEC, std::vector<std::vector<std::pair<int, int>>>>
		integrate_samples(const std::vector<STRVEC>& sss,
										const STRVEC& orig_samples);
	// 重複したサンプルが一つになるようにVCFを統合する
	static VCFSmall *integrate(const VCFFillable *vcf,
					const std::vector<std::vector<VCFFillableRecord *>>& rss,
					const STRVEC& orig_samples);
	static VCFSmall *merge(const std::vector<VCFFillable *>& vcfs,
							const STRVEC& orig_samples, const Option *option);
};
#endif
