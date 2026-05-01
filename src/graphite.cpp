#include <iostream>

#include "../include/graphite.h"
#include "../include/ImputeByRef.h"
#include "../include/ImputeNoRef.h"
#include "../include/materials.h"
#include "../include/option.h"
#include "../include/error_codes.h"


//////////////////// graphite ////////////////////

void graphite(const Option& option) {
	option.print_info();
	Materials	*materials = Materials::create(option.path_map,
													option.path_ped);
	materials->display_map_info();
	
	VCFHuge	*vcf = NULL;
	try {
		vcf = VCFHuge::read(option.path_vcf);
		if(option.exists_ref())
			ImputeByRef::impute(vcf, materials, option);
		else
			ImputeNoRef::impute(vcf, materials, option);
		delete vcf;
	}
	catch(const ExceptionWithCode& e) {
		if(vcf != NULL)
			delete vcf;
		delete materials;
		throw;
	}
	delete materials;
}

int main(int argc, char **argv) {
	const Option	*option = Option::create(argc, argv);
	if(option == NULL) {
		Option::usage(argv);
		exit(ErrorCode::WRONG_ARGUMENT);
	}
	
	try {
		graphite(*option);
	}
	catch(const ExceptionWithCode& e) {
		delete option;
		std::cerr << e.what() << std::endl;
		exit(e.get_error_code());
	}
	
	delete option;
	return 0;
}
