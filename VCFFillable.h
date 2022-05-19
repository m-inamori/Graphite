#ifndef __VCFFILLABLE
#define __VCFFILLABLE

#include "VCFFamily.h"


//////////////////// Genotype ////////////////////

class Genotype {
	const char	gt1;
	const char	gt2;
	const bool	phasing;
	
public:
	Genotype(const std::string& s);
	
	std::pair<char,char> gts() const;
	bool includes(char gt) const;
	bool is_consistent(const Genotype& mat_gt, const Genotype& pat_gt,
												bool considers_phasing) const;
	
	static bool is_valid(const std::string& gt);
	static int sum_gt(const std::string& gt);
};


//////////////////// VCFFillableRecord ////////////////////

class VCFFillableRecord : public VCFFamilyRecord {
public:
	enum class RecordType { UNABLE, FILLED, MAT, PAT, NONE };
	
protected:
	RecordType	type;
	
public:
	VCFFillableRecord(const STRVEC& v, const STRVEC& s, RecordType t) :
											VCFFamilyRecord(v, s), type(t) { }
	
	RecordType get_type() const { return type; }
	bool is_unable() const { return type == RecordType::UNABLE; }
	bool is_filled_type() const { return type == RecordType::FILLED; }
	bool is_mat_type() const { return type == RecordType::MAT; }
	bool is_pat_type() const { return type == RecordType::PAT; }
	
	STRVEC prog_gts() const;
	VCFFillableRecord *copy() const;
	std::string gt_from_parent(int mat_from, int pat_from) const;
	
	void modify();
	void modify_gts();
	int from_which_chrom(std::size_t i, bool is_mat) const;
	void modify_gts(const STRVEC& new_prog_gts);
	void fill_PGT();
	void set(const STRVEC& new_v, RecordType new_type);
	
private:
	std::vector<std::vector<double>> make_probability_table() const;
	int segregation_type() const;
	std::vector<std::size_t> conflicted_progeny_indices() const;
	// 例えば、子どもがだいたい0/0なのに親が0/0と1/1なら0/0に修正する
	// modifyできたかを返す
	bool modify_parents();
	bool is_valid();
	
	void disable() { this->type = RecordType::UNABLE; }
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
	
public:
	static VCFFillableRecord *from_VCFFamilyRecord(
									const VCFFamilyRecord *record);
	static VCFRecord *integrate_records(
						const std::vector<VCFFillableRecord *>& records);
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
		
		bool is_mat_prev_near() const;
		bool is_pat_prev_near() const;
		
		Pair select_nearest_froms(const std::vector<Pair>& pairs,
													std::size_t i) const;
		Pair select_pair(const std::vector<Pair>& pairs,
								std::size_t i, bool selected=true) const;
	};
	
	struct ConfigThread {
		const std::vector<VCFFillable *>&	vcfs;
		const std::size_t	first;
		const int	num_thread;
		
		ConfigThread(const std::vector<VCFFillable *>& v, int f, int n) :
								vcfs(v), first(f), num_thread(n) { }
	};
	
	using Position = std::tuple<int,ll,std::string>;
	using Group = std::pair<VCFFillableRecord::RecordType,
										std::vector<VCFFillableRecord *>>;
	
	std::vector<VCFFillableRecord *>	fillable_records;
	
public:
	VCFFillable(const std::vector<STRVEC>& h, const STRVEC& s,
									std::vector<VCFFillableRecord *> rs);
	~VCFFillable() { }
	
	std::vector<VCFRecord *> to_VCFRecord(std::vector<VCFFillableRecord *>& rs);
	VCFFillableRecord *get_record(std::size_t i) const {
		return fillable_records[i];
	}
	
	VCFFillable *insert_positions(const std::vector<Position>& positions);
	void impute();
	
private:
	std::vector<Group> group_records() const;
	VCFFillableRecord *find_prev_record(VCFFillableRecord::RecordType type,
							int i, const std::vector<Group>& groups) const;
	VCFFillableRecord *find_next_record(VCFFillableRecord::RecordType type,
							int i, const std::vector<Group>& groups) const;
	int from_which_chrom(const VCFFillableRecord *record,
										std::size_t i, bool is_mat) const;
	int from_which_chrom_mat(const VCFFillableRecord *record,
											std::size_t i) const {
		return from_which_chrom(record, i, true);
	}
	int from_which_chrom_pat(const VCFFillableRecord *record,
											std::size_t i) const {
		return from_which_chrom(record, i, false);
	}
	std::vector<double> probs_from_which_chrom(
								int prev_chrom, int next_chrom) const;
	double likelihood_each(const std::string& gt,
								const std::vector<double>& probs_mat,
								const std::vector<double>& probs_pat,
								int mat_phasing, int pat_phasing) const;
	double compute_phasing_likelihood_each(RecordSet& rs, std::size_t i,
									int mat_phasing, int pat_phasing) const;
	double compute_phasing_likelihood(RecordSet& rs, int mat_phasing,
														int pat_phasing) const;
	void phase(int i, const std::vector<Group>& groups);
	void determine_phasing(RecordSet& record_set);
	std::string modify_gt(RecordSet& rs, size_t i);
	void impute_core(RecordSet& record_set);
	
	VCFFillableRecord *find_prev_same_type_record(std::size_t i,
													std::size_t c) const;
	VCFFillableRecord *find_next_same_type_record(std::size_t i,
													std::size_t c) const;
	
	const RecordSet *create_recordset(std::size_t i,
										std::size_t c, bool is_mat) const;
	int select_mat(const std::vector<Pair>& pairs, const RecordSet *rs) const;
	void impute_NA_mat_each(std::size_t i, std::size_t c);
	void impute_NA_mat(std::size_t i);
	
	int select_pat(const std::vector<Pair>& pairs, const RecordSet *rs) const;
	void impute_NA_pat_each(std::size_t i, std::size_t c);
	void impute_NA_pat(std::size_t i);
	
public:
	static VCFFillable *convert(const VCFFamily *vcf);
	static void replace_filled_records(const std::vector<VCFFillable *>& vcfs,
															VCFHuge *orig_vcf);
	static std::vector<PosWithChr> merge_positions(
							std::vector<VCFFillable *>::const_iterator first,
							std::vector<VCFFillable *>::const_iterator last);
	static std::vector<PosWithChr> merge_positions(
									const std::vector<VCFFillable *>& vcfs);
	static STRVEC join_samples(const std::vector<VCFFillable *>& vcfs);
	static std::vector<STRVEC> join_header(
									const std::vector<VCFFillable *>& vcfs,
									const STRVEC& samples);
	static VCFSmall *merge_vcfs(const std::vector<VCFFillable *>& vcfs);
	static void impute_in_thread(void *config);
	static void impute_all_in_multithreads(
							const std::vector<VCFFillable *>& vcfs, int T);
	static VCFSmall *join_vcfs(const std::vector<VCFFamily *>& vcfs,
										const std::string& path_VCF, int T);
};
#endif
