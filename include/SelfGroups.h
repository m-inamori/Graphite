// SelfGroups.h
#ifndef __SELFGROUPS
#define __SELFGROUPS

#include <vector>

#include "VCFImpSelfRecord.h"
#include "SelfRecordSet.h"

class VCFSelfFillableRecord;


//////////////////// SelfGroups ////////////////////

class SelfGroups {
public:
	using Group = std::pair<SelfFillType, std::vector<VCFSelfFillableRecord *>>;
	
private:
	const std::vector<Group>	groups;
	
public:
	explicit SelfGroups(const std::vector<Group>& g) : groups(g) { }
	~SelfGroups() { }
	
	std::size_t size() const { return groups.size(); }
	const Group& get_group(std::size_t i) const { return groups[i]; }
	const SelfFillType get_type(std::size_t i) const { return groups[i].first; }
	const std::vector<VCFSelfFillableRecord *>&
	get_records(std::size_t i) const { return groups[i].second; }
	VCFSelfFillableRecord *find_prev_record(std::size_t i,
											SelfFillType g) const;
	VCFSelfFillableRecord *find_next_record(std::size_t i,
											SelfFillType g) const;
	std::vector<SelfRecordSet *> create_record_sets() const;
	
	static SelfGroups *create(const std::vector<VCFSelfFillableRecord *>& records);
};
#endif
