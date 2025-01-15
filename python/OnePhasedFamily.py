from __future__ import annotations

# coding: utf-8
# OnePhasedFamily.py
# 片親がphasingされて片親は分っているがphasingされていない家系を補完する

from functools import reduce
from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCFFamily import *
from VCFSmallFillable import *
from group import Groups
from RecordSet import RecordSet
from VCFImpFamily import FillType
from VCFHeteroHomoPP import *
from VCFImpHeteroHomo import *
from VCFHeteroImpHomo import *
from Map import *
import ClassifyRecord
from TypeDeterminer import *
from VCFOneParentImputed import VCFOneParentImputed
from VCFOneParentImputedRough import VCFOneParentImputedRough
from Genotype import Genotype


# VCFHeteroHomoPPを使わずにこれを使う あとで統合する
def classify_record(record: VCFRecord) -> tuple[ParentComb, FillType]:
	if record.is_NA(0) or record.is_NA(1):
		return (ParentComb.PNA, FillType.IMPUTABLE)
	
	i = 0 if record.v[9][0] != record.v[9][2] else 1
	j = 0 if record.v[10][0] != record.v[10][2] else 1
	if (i, j) == (0, 0):		# 0/1 x 0/1
		return (ParentComb.P01x01, FillType.IMPUTABLE)
	elif (i, j) == (1, 0):		# 0/0 x 0/1 or 1/1 x 0/1
		if record.v[9][0] == '0':
			return (ParentComb.P00x01, FillType.PAT)
		else:
			return (ParentComb.P01x11, FillType.PAT)
	elif (i, j) == (0, 1):		# 0/1 x 0/0 or 0/1 x 1/1
		if record.v[10][0] == '0':
			return (ParentComb.P00x01, FillType.MAT)
		else:
			return (ParentComb.P01x11, FillType.MAT)
	else:
		if record.v[9][0] == '0' and record.v[10][0] == '0':	# 0/0 x 0/0
			return (ParentComb.P00x00, FillType.FILLED)
		elif record.v[9][0] == '1' and record.v[10][0] == '1':	# 1/1 x 1/1
			return (ParentComb.P11x11, FillType.FILLED)
		else:													# 0/0 x 1/1
			return (ParentComb.P00x11, FillType.FILLED)

def classify_records(records: Sequence[VCFFamilyRecord]
								) -> list[list[VCFFillableRecord]]:
	# ヘテロ×ヘテロ, ホモ×ヘテロ, ヘテロ×ホモ, ホモ×ホモ
	rss: list[list[VCFFillableRecord]] = [ [] for _ in range(4) ]
	for index, record in enumerate(records):
		pair, type = classify_record(record)
		rss[type.value].append(VCFFillableRecord(record.v, record.samples,
														index, type, pair))
	return rss

def create(header: list[list[str]], rs: list[VCFFillableRecord],
						is_mat_hetero: bool, is_mat_imputed: bool,
						gmap: Map) -> VCFHeteroHomoOnePhased:
	if is_mat_hetero == is_mat_imputed:
		return VCFImpHeteroHomo(header, rs, is_mat_hetero, gmap)
	else:
		return VCFHeteroImpHomo(header, rs, is_mat_hetero, gmap)

def merge_vcf(rss: list[list[VCFFillableRecord]],
							header: list[list[str]]) -> VCFSmallFillable:
	rs = rss[0] + rss[1] + rss[2] + rss[3]
	rs.sort(key=lambda r: r.pos())
	return VCFSmallFillable(header, rs)

# HMMにrefを使っても計算量が十分に小さいか
def is_small(family: Family, ref_haps: list[list[int]]) -> bool:
	N = family.num_progenies()
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = NH**2 * 4**N * (2*NH + 2*N - 1)
	return R * M < 10**8 and R < 10**5

def is_small_ref(ref_haps: list[list[int]]) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = NH**2 * (2*NH - 1)
	return R * M < 10**8 and R < 10**5

def impute(family: Family, vcf: VCFFamily,
					unimputed_parents: list[str], gmap: Map) -> VCFSmallBase:
	header = vcf.header
	rss = classify_records(vcf.records)
	is_mat_imp = family.pat in unimputed_parents
	mat_vcf = create(header, rss[FillType.MAT.value], True, is_mat_imp, gmap)
	mat_vcf.impute()
	pat_vcf = create(header, rss[FillType.PAT.value], False, is_mat_imp, gmap)
	pat_vcf.impute()
	merged_vcf = merge_vcf(rss, header)
	merged_vcf.modify(False)
	return merged_vcf

def impute_by_parent(orig_vcf: VCFSmall, imputed_vcf: VCFSmall,
										ref_haps: list[list[int]],
										families: list[Family],
										non_imputed_parents: list[str],
										gmap: Map) -> VCFSmallBase:
	vcfs: list[VCFSmallBase] = []
	for family in families:
		vcf1 = VCFFamily.create_by_two_vcfs(imputed_vcf,
											orig_vcf, family.samples())
		is_mat_imputed = family.pat in non_imputed_parents
		if is_small(family, ref_haps):
			vcf = VCFOneParentImputed(vcf1.header, vcf1.records,
										ref_haps, is_mat_imputed, gmap)
			vcf.impute()
			vcfs.append(vcf)
		elif is_small_ref(ref_haps):
			vcf2 = VCFOneParentImputedRough(vcf1.header, vcf1.records,
										ref_haps, is_mat_imputed, gmap)
			vcf2.impute()
			vcfs.append(vcf2)
		else:
			imputed_vcf1 = impute(family, vcf1, non_imputed_parents, gmap)
			vcfs.append(imputed_vcf1)
	
	new_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute_by_parent']
