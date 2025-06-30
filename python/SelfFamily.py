from __future__ import annotations

# coding: utf-8
# SelfFamily.py
# 自殖の家系

from typing import Optional

from VCFFamily import *
from Map import *
from KnownFamily import *
from VCFSelfParentImputed import VCFSelfParentImputed
from VCFSelfProgenyImputed import VCFSelfProgenyImputed
from OptionSmall import OptionSmall


def impute(orig_vcf: VCFSmall, imputed_vcf: VCFSmall, ref_haps: list[list[int]],
									families: list[KnownFamily],
									imputed_samples: list[str],
									op: OptionSmall) -> Optional[VCFSmallBase]:
	vcfs: list[VCFSmallBase] = []
	set_imputed_samples = set(imputed_samples)
	for family in families:
		# imputedな後代のindex
		prog_indices = [ i for i, s in enumerate(family.progenies)
											if s in set_imputed_samples ]
		if family.mat in set_imputed_samples:
			# 親がimputedならimputedな後代はVCFに含めない
			samples = [family.mat] + [ s for s in family.progenies
											if s not in set_imputed_samples ]
			vcf = VCFSmall.create_by_two_vcfs(imputed_vcf, orig_vcf, samples)
			vcf1 = VCFSelfParentImputed(vcf.header, vcf.records, op.map)
			vcf1.impute()
			vcfs.append(vcf1)
		elif prog_indices:
			# 親がimputedでないならimputedな後代を一つだけVCFに含める
			samples = ([family.mat, family.progenies[prog_indices[0]]] +
									[ s for s in family.progenies
											if s not in set_imputed_samples ])
			vcf = VCFSmall.create_by_two_vcfs(imputed_vcf, orig_vcf, samples)
			vcf2 = VCFSelfProgenyImputed(vcf.header, vcf.records,
													ref_haps, 0, op.map)
			vcf2.impute()
			vcfs.append(vcf2)
	
	if not vcfs:
		return None
	
	print("%d self families have been imputed." % len(families))
	new_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
