from __future__ import annotations

# coding: utf-8
# ImputedAndKnownFamilyRef.py
# リファレンスがあって、
# 片親がphasingされて片親は分っているがphasingされていない家系を補完する

from typing import Optional, Sequence

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
import ImputedAndKnownFamily
import RefCommon
from OptionSmall import OptionSmall


#################### ImputedAndKnownFamilyRef ####################

def impute(orig_vcf: VCFSmall, phased_vcf: VCFGeno,
									ref_haps: list[list[int]],
									families: list[Family],
									non_imputed_parents: list[str],
									op: OptionSmall) -> Optional[VCFGenoBase]:
	vcfs: list[VCFGenoBase] = []
	for family in families:
		vcf = VCFGeno.extract_samples(family.samples(), orig_vcf);
		records = RefCommon.merge_family_records(phased_vcf,
													vcf, family.samples())
		is_mat_imputed = family.pat in non_imputed_parents
		vcf_family = ImputedAndKnownFamily.create_family_vcf(
										family, records, is_mat_imputed,
										len(families), ref_haps, orig_vcf, op)
		vcf_family.impute()
		vcfs.append(vcf_family)
	
	if not vcfs:
		return None
	
	print("%d families whose one parent is imputed and the other parent is"
									" known have been imputed." % len(vcfs))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
