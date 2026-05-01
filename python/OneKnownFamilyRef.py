from __future__ import annotations

# coding: utf-8
# OneImputedFamilyRef.py
# 片親がimputeされていてもう片親はunknownな家系を補完する

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
import OneKnownFamily
import RefCommon
from KnownFamily import *
from OptionSmall import OptionSmall


#################### OneImputedFamilyRef ####################

def impute(orig_vcf: VCFSmall, phased_vcf: VCFGeno, ref_haps: list[list[int]],
									families: list[KnownFamily],
									op: OptionSmall) -> Optional[VCFGenoBase]:
	vcfs: list[VCFGenoBase] = []
	for family in families:
		vcf = VCFGeno.extract_samples(family.samples(), orig_vcf);
		records = RefCommon.merge_family_records(phased_vcf,
													vcf, family.samples())
		vcf1 = OneKnownFamily.create_family_vcf(family, records,
													len(families),
													ref_haps, orig_vcf, op)
		if vcf1 is not None:
			vcf1.impute()
			vcfs.append(vcf1)
	
	if not vcfs:
		return None
	
	print("%d families whose one parent is known and the other parent is"
								" unknown have been imputed." % len(families))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
