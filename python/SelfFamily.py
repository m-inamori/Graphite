from __future__ import annotations

# coding: utf-8
# SelfFamily.py
# 自殖の家系

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFSelfImputable import VCFSelfImputable
from GenoRecord import GenoRecord
from Map import *
from KnownFamily import *
from VCFSelfParentImputed import VCFSelfParentImputed
from VCFSelfProgenyImputed import VCFSelfProgenyImputed
from OptionSmall import OptionSmall


#################### SelfFamily ####################

def create_record(record: GenoRecord, indices: list[int]) -> GenoRecord:
	pos = record.pos
	geno = [ record.geno[i] for i in indices ]
	return GenoRecord(pos, geno)

def create_family_vcf(orig_vcf: VCFSmall, records: list[GenoRecord],
								ref_haps: list[list[int]],
								family: KnownFamily,
								imputed_samples: list[str],
								op: OptionSmall) -> Optional[VCFSelfImputable]:
	set_imputed_samples = set(imputed_samples)
	# imputedな後代のindex
	prog_indices = [ i for i, s in enumerate(family.progenies)
										if s in set_imputed_samples ]
	if family.mat in set_imputed_samples:
		# 親がimputedならimputedな後代はVCFに含めない
		indices = [0] + [ i+1 for i, s in enumerate(family.progenies)
										if s not in set_imputed_samples ]
		samples = [family.mat] + [ s for s in family.progenies
										if s not in set_imputed_samples ]
		new_records: list[GenoRecord] = [ create_record(r, indices)
														for r in records ]
		return VCFSelfParentImputed(samples, new_records, op.map, orig_vcf)
	elif prog_indices:
		# 親がimputedでないならimputedな後代を一つだけVCFに含める
		samples = ([family.mat, family.progenies[prog_indices[0]]] +
								[ s for s in family.progenies
										if s not in set_imputed_samples ])
		indices = ([0, prog_indices[0]+1] + 
								[ i for i, s in enumerate(family.progenies)
											if s not in set_imputed_samples ])
		new_records = [ create_record(r, indices) for r in records ]
		return VCFSelfProgenyImputed(samples, new_records,
											ref_haps, 0, op.map, orig_vcf)
	else:
		return None

def impute(orig_vcf: VCFSmall, imputed_vcf: VCFGeno, ref_haps: list[list[int]],
										families: list[KnownFamily],
										imputed_samples: list[str],
										op: OptionSmall) -> Optional[VCFGeno]:
	vcfs: list[VCFGenoBase] = []
	set_imputed_samples = set(imputed_samples)
	for family in families:
		vcf = VCFGeno.create_by_two_vcfs(imputed_vcf,
											orig_vcf, family.samples())
		vcf_family = create_family_vcf(orig_vcf, vcf.records, ref_haps,
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
