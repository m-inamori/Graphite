from __future__ import annotations

# coding: utf-8
# OnePhasedFamily.py
# 片親がphasingされて片親は分っているがphasingされていない家系を補完する

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

def impute(family: Family, vcf: VCFFamily,
					non_imputed_parents: list[str], gmap: Map) -> VCFSmallBase:
	header = vcf.header
	is_mat_imp = family.pat in non_imputed_parents
	N = family.num_progenies()
	M = len(vcf)
	if 4**(N*2+1)*N*M <= 10**8:
		vcfopi = VCFOneParentImputed(header, vcf.records, is_mat_imp, gmap)
		vcfopi.impute()
		return
	
	rss = classify_records(vcf.records)
	mat_vcf = create(header, rss[FillType.MAT.value], True, is_mat_imp, gmap)
	mat_vcf.impute()
	pat_vcf = create(header, rss[FillType.PAT.value], False, is_mat_imp, gmap)
	pat_vcf.impute()
	merged_vcf = merge_vcf(rss, header)
	merged_vcf.modify(False)
	return merged_vcf

def impute_by_parent(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
									families: list[Family],
									unimputed_parents: list[str],
									gmap: Map) -> VCFSmallBase:
	vcfs: list[VCFSmallBase] = []
	for family in families:
		vcf = VCFFamily.create_by_two_vcfs(merged_vcf,
											orig_vcf, family.samples())
		impute(family, vcf, unimputed_parents, gmap)
		vcfs.append(vcf)
	
	new_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute_by_parent']
