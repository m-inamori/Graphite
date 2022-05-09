#ifndef __PEDIGREE
#define __PEDIGREE

#include <vector>
#include <string>
#include <set>
#include <utility>

class VCFBase;


//////////////////// Progeny ////////////////////

class Progeny {
	const std::string	family;
	const std::string	name;
	const std::string	mat;
	const std::string	pat;
	
public:
	Progeny(const std::string& f, const std::string& n,
			const std::string& m, const std::string& p) :
						family(f), name(n), mat(m), pat(p) { }
	
	const Progeny *copy() const;
	std::string get_name() const { return name; }
	std::pair<std::string, std::string> parents() const {
		return std::pair<std::string, std::string>(mat, pat);
	}
	bool with_parents() const;
	bool is_all_in_samples(const std::set<std::string>& samples) const;
	
public:
	static const Progeny *create(const std::vector<std::string>& v);
};


//////////////////// PedigreeTable ////////////////////

class PedigreeTable {
	const std::vector<const Progeny *>	table;
	
public:
	PedigreeTable(const std::vector<const Progeny *>& progs) : table(progs) { }
	~PedigreeTable();
	
	std::size_t size() const { return table.size(); }
	const PedigreeTable *filter_with_parents(
						const std::vector<std::string>& samples) const;
	std::vector<std::pair<std::string,std::string>> extract_families() const;
	std::vector<std::string> get_family_children(const std::string& mat,
												 const std::string& pat) const;
	
public:
	static const PedigreeTable *read(const std::string& path);
	static const PedigreeTable *create(const std::string& path, VCFBase *vcf);
};
#endif
