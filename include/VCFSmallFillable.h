#ifndef __VCFSMALLFILLABLE
#define __VCFSMALLFILLABLE

#include "VCFFillable.h"

class Map;


//////////////////// VCFSmallFillable ////////////////////

class VCFSmallFillable : public VCFFillable {
public:
	struct ConfigThreadPhase {
		const std::size_t	first;
		const std::size_t	num_threads;
		const std::vector<RecordSetSmall *>&	record_sets;
		VCFFillable	*vcf;
		
		ConfigThreadPhase(std::size_t f, std::size_t n,
					const std::vector<RecordSetSmall *>& r, VCFFillable *v) :
					first(f), num_threads(n), record_sets(r), vcf(v) { }
	};
	
public:
	VCFSmallFillable(const std::vector<STRVEC>& h, const STRVEC& s,
							const std::vector<VCFFillableRecord *>& rs) :
													VCFFillable(h, s, rs) { }
	~VCFSmallFillable() { }
	
	///// non-virtual methods /////
	void delete_records();
	void modify(int T);
	
	const RecordSet *create_recordset(
							std::size_t i, std::size_t c, bool is_mat) const;
	
public:
	static void phase_in_thread(void *config);
};
#endif
