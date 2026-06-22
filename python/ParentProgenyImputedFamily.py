from __future__ import annotations

# coding: utf-8
# ParentProgenyImputedFamily.py
# 片親と後代がphasingされて片親は分っているがphasingされていない家系を補完する

from typing import Optional, Sequence

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import VCFFamily, VCFFamilyRecord
from VCFImputable import *
from VCFOneParentProgenyImputed import VCFOneParentProgenyImputed
from KnownFamily import KnownFamily
from Map import *
from OptionSmall import OptionSmall


#################### ParentProgenyImputedFamily ####################

# HMMにrefを使っても計算量が十分に小さいか
def is_small(ref_haps: list[list[int]],
					L: int, P: int, op: OptionSmall) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = (NH * (NH+P*2+1) << (P*2+1)) / op.precision_ratio
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

# 適正な後代の個数を調べる
def compute_upper_num_progenies(ref_haps: list[list[int]],
								L: int, P: int, op: OptionSmall) -> int:
	for p in range(P, -1, -1):
		if is_small(ref_haps, L, p, op):
			return p
	else:
		return 0	# dummy

def reduce_progenies(p: int, f: KnownFamily) -> KnownFamily:
	progs = f.progenies[:p]
	return KnownFamily(f.mat, f.pat, f.mat_known, f.pat_known, progs)

def create_family_vcf(family: KnownFamily, records: list[VCFFamilyRecord],
						num_families: int, ref_haps: list[list[int]],
						should_impute_mat: bool, orig_vcf: VCFSmall,
						op: OptionSmall) -> VCFImputable:
	vcf = VCFFamily(family.samples(), records, orig_vcf)
	num_progs = len(family.samples()) - 2
	p = compute_upper_num_progenies(ref_haps, num_families, num_progs, op)
	if p < num_progs:
		family = reduce_progenies(p, family)
		vcf = vcf.extract(family.samples())
	
	return VCFOneParentProgenyImputed(family.samples(), vcf.records,
										ref_haps, should_impute_mat,
										op.map, orig_vcf)

def impute(orig_vcf: VCFSmall, imputed_vcf: VCFGeno,
									ref_haps: list[list[int]],
									families: list[KnownFamily],
									non_imputed_parents: list[str],
									op: OptionSmall) -> Optional[VCFGenoBase]:
	vcfs: list[VCFGenoBase] = []
	for family in families:
		vcf = VCFFamily.create_by_two_vcfs(imputed_vcf,
											orig_vcf, family.samples())
		should_impute_mat = vcf.mat() in non_imputed_parents
		parent = vcf.mat() if should_impute_mat else vcf.pat()
		vcf1 = create_family_vcf(family, vcf.records, len(families),
											ref_haps, should_impute_mat,
											orig_vcf, op)
		vcf1.impute()
		imputed_samples = [parent] + vcf1.samples[3:]
		vcf2 = vcf1.extract_by_samples(imputed_samples)
		vcfs.append(vcf2)
	
	if not vcfs:
		return None
	
	print("%d families whose one parent and"
					" one progeny is imputed have been imputed." % len(vcfs))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute', 'impute_family']
