from __future__ import annotations

# coding: utf-8
# graphite.py

from itertools import *
import sys

import ImputeByRef
import ImputeNoRef
from VCF import VCFHuge
from materials import *
from Map import *
from option import *
from exception_with_code import *


#################### process ####################

def graphite(option: Option) -> None:
	option.print_info()
	materials = Materials.create(option.path_map, option.path_ped)
	materials.display_map_info()
	vcf = VCFHuge.read(option.path_VCF)
	
	if option.exists_ref():
		ImputeByRef.impute(vcf, materials, option)
	else:
		ImputeNoRef.impute(vcf, materials, option)


#################### main ####################

option = Option.create(sys.argv)
try:
	if option is None:
		Option.usage()
		exit(1)
	else:
		graphite(option)
		exit(0)
except ExceptionWithCode as e:
	print(str(e), file=sys.stderr)
	exit(e.get_error_code().value)
