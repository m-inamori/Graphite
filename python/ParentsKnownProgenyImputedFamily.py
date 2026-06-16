from __future__ import annotations

# coding: utf-8
# ParentsKnownProgenyImputedFamily.py
# 両親がknownで後代が一つ以上imputeされている家系を補完する

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
from VCFProgenyImputed import *
from VCFParentsKnownProgenyImputed import *
from Map import *
from KnownFamily import KnownFamily
from ReferenceHaplotype import filter_haplotypes
from OptionSmall import OptionSmall


#################### ParentsKnownProgenyImputedFamily ####################

def is_small(ref_haps: list[list[int]], L: int, op: OptionSmall) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	# 4つの一つだけ変わってもOK
	R = NH**2 * (8*NH + 4) / op.precision_ratio
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

# is_smallが通るギリギリのNHを求める
def compute_upper_NH(family: Family, M: int, L: int, op: OptionSmall) -> int:
	for NH in count(2):
		R = NH**2 * (8*NH + 4) / op.precision_ratio
		if not (R * M < 10**8 and R < 10**5 and L * R * M < 10**9):
			break
	return NH - 1

def impute_family(family: KnownFamily, samples: list[str],
					records: list[VCFFamilyRecord],
					num_families: int, ref_haps: list[list[int]],
					vcf: VCFSmall, op: OptionSmall) -> VCFGenoBase:
	NH = compute_upper_NH(family, len(records), num_families, op)
	if is_small(ref_haps, num_families, op):
		vcf1 = VCFParentsKnownProgenyImputed(samples, records,
												ref_haps, ref_haps,
												op.map, vcf)
		vcf1.impute()
		return vcf1
	elif NH >= 10:
		NH2 = min(20, NH)
		gts_mat = [ r.geno[0] for r in records ]
		gts_pat = [ r.geno[1] for r in records ]
		ref_haps_mat = filter_haplotypes(ref_haps, gts_mat, NH2)
		ref_haps_pat = filter_haplotypes(ref_haps, gts_pat, NH2)
		vcf2 = VCFParentsKnownProgenyImputed(samples, records,
												ref_haps_mat, ref_haps_pat,
												op.map, vcf)
		vcf2.impute()
		return vcf2
	else:
		# 計算量が大きいときは、片親だけ補完する
		new_samples = [samples[0], samples[2]]
		new_records = [ VCFFamilyRecord(r.pos, [r.geno[0], r.geno[2]])
														for r in records ]
		vcf3 = VCFProgenyImputed(new_samples, new_records, ref_haps,
													True, op.map, vcf)
		vcf3.impute()
		# 補完した親だけのVCFにして返す
		vcf_parent = VCFGenoBase.extract_by_samples(vcf3, samples[:1])
		return vcf_parent

def impute(orig_vcf: VCFSmall, merged_vcf: VCFGenoBase,
			families: list[KnownFamily], imputed_progenies: list[list[str]],
			ref_haps: list[list[int]], op: OptionSmall) -> Optional[VCFGeno]:
	vcfs: list[VCFGenoBase] = []
	for i in range(len(families)):
		family = families[i]
		parents = [family.mat, family.pat]
		progeny = imputed_progenies[i][0]
		samples = parents + [progeny]
		vcf = VCFFamily.create_by_two_vcfs(merged_vcf, orig_vcf, samples)
		vcf1 = impute_family(family, samples, vcf.records,
									len(families), ref_haps, orig_vcf, op)
		vcfs.append(vcf1)
	
	if not vcfs:
		return None
	
	print("%d families whose progeny is imputed have been imputed." %
															len(families))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
