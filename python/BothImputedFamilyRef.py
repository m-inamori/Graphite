from __future__ import annotations

# coding: utf-8
# BothImputedFamilyRef.py
# リファレンスがあって、両親がimputeされている家系を補完する

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
from VCFBothParentImputed import *
from Map import *
from KnownFamily import KnownFamily
import RefCommon
from OptionSmall import OptionSmall


#################### BothImputedFamilyRef ####################

def impute(orig_vcf: VCFSmall, phased_vcf: VCFGenoBase,
			families: list[KnownFamily], op: OptionSmall) -> Optional[VCFGeno]:
	vcfs: list[VCFGenoBase] = []
	for family in families:
		vcf = VCFGeno.extract_samples(family.samples(), orig_vcf);
		records = RefCommon.merge_family_records(phased_vcf,
													vcf, family.samples())
		family_vcf = VCFBothParentImputed(vcf.samples, records,
														op.map, vcf.vcf)
		family_vcf.impute(op.precision_ratio, op.num_threads)
		vcfs.append(family_vcf)
	
	if not vcfs:
		return None
	
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	print("%d both parent imputed families have been imputed." % len(vcfs))
	return new_vcf

__all__ = ['impute']
