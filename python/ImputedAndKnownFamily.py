from __future__ import annotations

# coding: utf-8
# ImputedAndKnownFamily.py
# 片親がphasingされて片親は分っているがphasingされていない家系を補完する

from typing import Optional, Sequence

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
from VCFSmallFillable import *
from group import Groups
from RecordSet import RecordSet
from VCFImpFamilyRecord import FillType
from VCFHeteroHomoPP import *
from VCFImpHeteroHomo import *
from VCFHeteroImpHomo import *
from Map import *
import ClassifyRecord
from TypeDeterminer import *
from VCFOneParentImputed import VCFOneParentImputed
from VCFOneParentImputedRough import VCFOneParentImputedRough
from OptionSmall import OptionSmall


# VCFHeteroHomoPPを使わずにこれを使う あとで統合する
def classify_record(record: VCFFamilyRecord) -> tuple[ParentComb, FillType]:
	if record.is_mat_NA() or record.is_pat_NA():
		return (ParentComb.PNA, FillType.IMPUTABLE)
	
	i = 0 if record.is_mat_hetero() else 1
	j = 0 if record.is_pat_hetero() else 1
	if (i, j) == (0, 0):		# 0/1 x 0/1
		return (ParentComb.P01x01, FillType.IMPUTABLE)
	elif (i, j) == (1, 0):		# 0/0 x 0/1 or 1/1 x 0/1
		if record.is_00(0):
			return (ParentComb.P00x01, FillType.PAT)
		else:
			return (ParentComb.P01x11, FillType.PAT)
	elif (i, j) == (0, 1):		# 0/1 x 0/0 or 0/1 x 1/1
		if record.is_00(1):
			return (ParentComb.P00x01, FillType.MAT)
		else:
			return (ParentComb.P01x11, FillType.MAT)
	else:
		if record.is_00(0) and record.is_00(1):
			return (ParentComb.P00x00, FillType.FILLED)
		elif record.is_11(0) and record.is_11(1):
			return (ParentComb.P11x11, FillType.FILLED)
		else:													# 0/0 x 1/1
			return (ParentComb.P00x11, FillType.FILLED)

def classify_records(samples: list[str],
						records: Sequence[VCFFamilyRecord],
						ref_vcf: VCFSmall) -> list[list[VCFFillableRecord]]:
	cols = ref_vcf.extract_columns(samples)
	# ヘテロ×ヘテロ, ホモ×ヘテロ, ヘテロ×ホモ, ホモ×ホモ
	rss: list[list[VCFFillableRecord]] = [ [] for _ in range(4) ]
	for index, record in enumerate(records):
		pair, type = classify_record(record)
		ref_record = ref_vcf.records[index]
		probs = ref_record.parse_PL(record.geno, cols)
		rss[type.value].append(VCFFillableRecord(record.pos, record.geno,
													index, type, pair, probs))
	return rss

def create(samples: list[str], rs: list[VCFFillableRecord],
						is_mat_hetero: bool, is_mat_imputed: bool,
						gmap: Map, vcf: VCFSmall) -> VCFHeteroHomoOnePhased:
	if is_mat_hetero == is_mat_imputed:
		return VCFImpHeteroHomo(samples, rs, is_mat_hetero, gmap, vcf)
	else:
		return VCFHeteroImpHomo(samples, rs, is_mat_hetero, gmap, vcf)

def merge_vcf(samples: list[str], rss: list[list[VCFFillableRecord]],
												vcf: VCFSmall) -> VCFFillable:
	rs = rss[0] + rss[1] + rss[2] + rss[3]
	rs.sort(key=lambda r: r.pos)
	return VCFSmallFillable(samples, rs, vcf)

# HMMにrefを使っても計算量が十分に小さいか
def is_small(family: Family, ref_haps: list[list[int]],
									L: int, op: OptionSmall) -> bool:
	N = family.num_progenies()
	if N > 2:
		return False
	
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R: float = NH**2 * 4**N * (2*NH + 2*N - 1) / op.precision_ratio
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

def is_small_ref(ref_haps: list[list[int]], L: int, op: OptionSmall) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = NH**2 * (2*NH - 1) / op.precision_ratio
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

def impute_roughly(family: Family, vcf: VCFFamily,
					unimputed_parents: list[str], gmap: Map) -> VCFFillable:
	samples = family.samples()
	rss = classify_records(samples, vcf.records, vcf.vcf)
	is_mat_imp = family.pat in unimputed_parents
	mat_vcf = create(samples, rss[FillType.MAT.value], True,
											is_mat_imp, gmap, vcf.vcf)
	mat_vcf.impute()
	pat_vcf = create(samples, rss[FillType.PAT.value], False,
											is_mat_imp, gmap, vcf.vcf)
	pat_vcf.impute()
	merged_vcf = merge_vcf(samples, rss, vcf.vcf)
	merged_vcf.modify(False)
	return merged_vcf

def impute(orig_vcf: VCFSmall, imputed_vcf: VCFGeno,
									ref_haps: list[list[int]],
									families: list[Family],
									non_imputed_parents: list[str],
									op: OptionSmall) -> Optional[VCFGenoBase]:
	vcfs: list[VCFGenoBase] = []
	for family in families:
		vcf1 = VCFFamily.create_by_two_vcfs(imputed_vcf,
											orig_vcf, family.samples())
		is_mat_imputed = family.pat in non_imputed_parents
		if is_small(family, ref_haps, len(families), op):
			vcf = VCFOneParentImputed(vcf1.samples, vcf1.records, ref_haps,
											is_mat_imputed, op.map, vcf1.vcf)
			vcf.impute()
			vcfs.append(vcf)
		elif is_small_ref(ref_haps, len(families), op):
			vcf2 = VCFOneParentImputedRough(vcf1.samples, vcf1.records,
									ref_haps, is_mat_imputed, op.map, vcf1.vcf)
			vcf2.impute()
			vcfs.append(vcf2)
		else:
			imputed_vcf1 = impute_roughly(family, vcf1,
											non_imputed_parents, op.map)
			vcfs.append(imputed_vcf1)
	
	if not vcfs:
		return None
	
	print("%d families whose one parent is imputed and the other parent is"
									" known have been imputed." % len(vcfs))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
