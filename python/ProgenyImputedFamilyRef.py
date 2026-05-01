from __future__ import annotations

# coding: utf-8
# ProgenyImputedFamilyRef.py
# 片親がknownで後代が一つ以上imputeされている家系を補完する

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
import ProgenyImputedFamily
import RefCommon
from KnownFamily import *
from OptionSmall import OptionSmall


#################### ProgenyImputedFamilyRef ####################

def impute(orig_vcf: VCFSmall, phased_vcf: VCFGenoBase,
			families: list[KnownFamily], imputed_progenies: list[list[str]],
			ref_haps: list[list[int]], op: OptionSmall) -> Optional[VCFGeno]:
	vcfs: list[VCFGenoBase] = []
	for i in range(len(families)):
		family = families[i]
		parent = family.mat if family.mat_known else family.pat
		progeny = imputed_progenies[i][0]
		samples = [parent, progeny]
		vcf = VCFGeno.extract_samples(samples, orig_vcf);
		records = RefCommon.merge_family_records(phased_vcf, vcf, samples)
		vcf1 = ProgenyImputedFamily.impute_family(family, samples, records,
													len(families), ref_haps,
													orig_vcf, op)
		if vcf1 is not None:
			vcfs.append(vcf1)
	
	if not vcfs:
		return None
	
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	
	print("%d self families have been imputed." % len(families))
	return new_vcf

__all__ = ['impute']
