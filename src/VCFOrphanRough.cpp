#include <cmath>
#include "../include/VCFOrphanRough.h"
#include "../include/OrphanImputer.h"
#include "../include/common.h"

using namespace std;

VCFOrphanRough::VCFOrphanRough(
				const STRVEC& s, const std::vector<GenoRecord *>& rs,
				const std::vector<std::vector<std::vector<int>>>& ref_hs_table,
				const Map& map_, double w, const VCFSmall *vcf) :
							VCFGenoBase(s, vcf),
							records(rs),
							imputers(create_imputers(ref_hs_table, map_, w)) { }

VCFOrphanRough::~VCFOrphanRough() {
	Common::delete_all(records);
	Common::delete_all(imputers);
}

vector<OrphanImputer *> VCFOrphanRough::create_imputers(
							const vector<vector<vector<int>>>& ref_haps_table,
							const Map& map_, double w) {
	vector<OrphanImputer *>	imputers;
	for(auto p = ref_haps_table.begin(); p != ref_haps_table.end(); ++p) {
		imputers.push_back(new OrphanImputer(records, *p, map_, w));
	}
	return imputers;
}

void VCFOrphanRough::impute(size_t io) {
	imputers[io]->impute(io);
}
