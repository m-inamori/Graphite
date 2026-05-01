from __future__ import annotations

# coding: utf-8
# ProgenyImputedFamily.py
# 片親がknownで後代が一つ以上imputeされている家系を補完する

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
from VCFProgenyImputed import *
from Map import *
from KnownFamily import KnownFamily
from OptionSmall import OptionSmall


#################### ProgenyImputedFamily ####################

def impute_family(family: KnownFamily, samples: list[str],
					records: list[VCFFamilyRecord],
					num_families: int, ref_haps: list[list[int]],
					vcf: VCFSmall, op: OptionSmall) -> VCFGenoBase:
	vcf1 = VCFProgenyImputed(samples, records, ref_haps,
										family.mat_known, op.map, vcf)
	vcf1.impute()
	return vcf1

def impute(orig_vcf: VCFSmall, merged_vcf: VCFGenoBase,
			families: list[KnownFamily], imputed_progenies: list[list[str]],
			ref_haps: list[list[int]], op: OptionSmall) -> Optional[VCFGeno]:
	vcfs: list[VCFGenoBase] = []
	for i in range(len(families)):
		family = families[i]
		parent = family.mat if family.mat_known else family.pat
		progeny = imputed_progenies[i][0]
		samples = [parent, progeny]
		vcf = VCFFamily.create_by_two_vcfs(merged_vcf, orig_vcf, samples)
		vcf1 = impute_family(family, samples, vcf.records,
									len(families), ref_haps, orig_vcf, op)
		vcfs.append(vcf1)
	
	if not vcfs:
		return None
	
	print("%d families whose progeny is imputed have been imputed." %
															len(families))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
