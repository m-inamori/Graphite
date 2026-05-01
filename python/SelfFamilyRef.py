from __future__ import annotations

# coding: utf-8
# SelfFamily.py
# 自殖の家系

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
import SelfFamily
import RefCommon
from KnownFamily import *
from OptionSmall import OptionSmall


def impute(orig_vcf: VCFSmall, phased_vcf: VCFGeno, ref_haps: list[list[int]],
										families: list[KnownFamily],
										imputed_samples: list[str],
										op: OptionSmall) -> Optional[VCFGeno]:
	vcfs: list[VCFGenoBase] = []
	set_imputed_samples = set(imputed_samples)
	for family in families:
		vcf = VCFGeno.extract_samples(family.samples(), orig_vcf);
		records = RefCommon.merge_records(phased_vcf, vcf, family.samples())
		vcf_family = SelfFamily.create_family_vcf(orig_vcf, records, ref_haps,
													family, imputed_samples, op)
		if vcf_family is not None:
			vcf_family.impute()
			vcfs.append(vcf_family)
	
	if not vcfs:
		return None
	
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	
	print("%d self families have been imputed." % len(families))
	return new_vcf

__all__ = ['impute']
