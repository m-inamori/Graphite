from __future__ import annotations

# coding: utf-8
# ProgenyImputedFamilyRef.py
# 片親がknownで後代が一つ以上imputeされている家系を補完する

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFProgenyImputed import VCFProgenyImputed
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
		progeny = imputed_progenies[i][0]
		samples = [family.mat, family.pat, progeny]
		vcf = VCFGeno.extract_samples(samples, orig_vcf);
		records = RefCommon.merge_family_records(phased_vcf, vcf, samples)
		vcf1 = VCFProgenyImputed(samples, records, ref_haps,
									family.mat_known, op.map, orig_vcf)
		vcf1.impute()
		parent = family.mat if family.mat_known else family.pat
		vcf_parent = VCFGenoBase.extract_by_samples(vcf1, [parent])
		vcfs.append(vcf_parent)
	
	if not vcfs:
		return None
	
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	
	print("%d self families have been imputed." % len(families))
	return new_vcf

__all__ = ['impute']
