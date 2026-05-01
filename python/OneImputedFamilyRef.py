from __future__ import annotations

# coding: utf-8
# OneImputedFamilyRef.py
# 片親がimputeされていてもう片親はunknownな家系を補完する

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
import OneImputedFamily
import RefCommon
from KnownFamily import *
from OptionSmall import OptionSmall


#################### OneImputedFamilyRef ####################

def impute(orig_vcf: VCFSmall, phased_vcf: VCFGeno, ref_haps: list[list[int]],
									families: list[KnownFamily],
									op: OptionSmall) -> Optional[VCFGenoBase]:
	vcfs: list[VCFGenoBase] = []
	set_imputed_samples = set(phased_vcf.samples)
	for family in families:
		is_mat_imputed = family.mat in set_imputed_samples
		vcf = VCFGeno.extract_samples(family.samples(), orig_vcf);
		records = RefCommon.merge_family_records(phased_vcf,
													vcf, family.samples())
		vcf1 = OneImputedFamily.create_family_vcf(family, is_mat_imputed,
													records, len(families),
													ref_haps, orig_vcf, op)
		vcf1.impute()
		vcfs.append(vcf1)
	
	if not vcfs:
		return None
	
	print("%d families whose one parent is imputed and the other parent is"
								" unknown have been imputed." % len(families))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
