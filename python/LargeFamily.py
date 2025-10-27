from __future__ import annotations

# coding: utf-8
# LargeFamily.py

from itertools import *
from collections import defaultdict, Counter
import sys
from multiprocessing import Pool

from typing import Optional

from VCF import *
from VCFGeno import VCFGeno
from VCFFamily import VCFFamily, VCFFamilyRecord
from pedigree import PedigreeTable, Family
from VCFHeteroHomo import *
from VCFHomoHomo import VCFHomoHomoRecord
from VCFHeteroHeteroLiteRecord import VCFHeteroHeteroLiteRecord
from VCFJunkRecord import *
from VCFFillable import *
from VCFHeteroHomoPP import *
from VCFImpFamilyRecord import VCFImpFamilyRecord
from Map import *
import ClassifyRecord as CR
from KnownFamily import KnownFamily
from TypeDeterminer import TypeDeterminer, ParentComb
from Genotype import Genotype
from option import *
from common import *


#################### process ####################

def create_heterohomo_record(pos: int, geno: list[int], family: KnownFamily,
			i: int, wrong_type: str, pair: ParentComb) -> VCFHeteroHomoRecord:
	# 片親が不明の時、Genotypeを補う
	total_int_gt = 1 if pair == ParentComb.P00x01 else 3
	if not family.mat_known:
		pat_int_gt = Genotype.unphased(geno[1])
		mat_int_gt = total_int_gt - pat_int_gt
		if pat_int_gt == 3 or mat_int_gt < 0 or 2 < mat_int_gt:
			wrong_type = 'Modifiable'
		else:
			geno[0] = mat_int_gt
	elif not family.pat_known:
		mat_int_gt = Genotype.unphased(geno[0])
		pat_int_gt = total_int_gt - mat_int_gt
		if mat_int_gt == 3 or pat_int_gt < 0 or 2 < pat_int_gt:
			wrong_type = 'Modifiable'
		else:
			geno[1] = pat_int_gt
	return VCFHeteroHomoRecord(pos, geno, i, wrong_type, pair)

def classify_record(i: int, vcf: VCFFamily, td: TypeDeterminer,
					family: KnownFamily,
					heho_records: list[Optional[VCFHeteroHomoRecord]],
					other_records: list[Optional[VCFImpFamilyRecord]]) -> None:
	record = vcf.records[i]
	pos = record.pos
	geno = vcf.records[i].geno
	samples = vcf.samples
	wrong_type, pair = CR.classify_record(record, td, family.is_one_unknown())
	if pair.is_homohomo():
		record_ = VCFHomoHomoRecord(pos, geno, i, wrong_type, pair)
		other_records[i] = record_.impute()
	elif pair.is_heterohomo():
		heho_records[i] = create_heterohomo_record(pos, geno, family, i,
														wrong_type, pair)
	elif pair == ParentComb.P01x01:
		other_records[i] = VCFHeteroHeteroLiteRecord(pos, geno, i, wrong_type)
	else:		# 候補が無い
		other_records[i] = VCFJunkRecord(pos, geno, i, wrong_type)

def classify_records_parallel(v: tuple[VCFFamily, TypeDeterminer, KnownFamily,
									   list[Optional[VCFHeteroHomoRecord]],
									   list[Optional[VCFImpFamilyRecord]],
									   int, int]) -> None:
	vcf, td, family, heho_records, other_records, i0, num_threads = v
	for i in range(i0, len(vcf.records), num_threads):
		classify_record(i, vcf, td, family, heho_records, other_records)

# HeteroHomoだけ別にする
# このあとHeteroHomoだけ補完するから
# その他はVCFFillableにした後補完する
def classify_records(vcf: VCFFamily, family: KnownFamily,
						option: Option) -> tuple[list[VCFHeteroHomoRecord],
												 list[VCFImpFamilyRecord]]:
	td = CR.get_typedeterminer(vcf.num_progenies(), option.ratio)
	records1: list[Optional[VCFHeteroHomoRecord]] = [None] * len(vcf)
	records2: list[Optional[VCFImpFamilyRecord]] = [None] * len(vcf)
	args = [ (vcf, td, family, records1, records2, i, option.num_threads)
										for i in range(option.num_threads) ]
	for arg in args:
		classify_records_parallel(arg)
	
	heho_records: list[VCFHeteroHomoRecord] = []
	other_records: list[VCFImpFamilyRecord] = []
	for i in range(len(vcf)):
		r1 = records1[i]
		r2 = records2[i]
		if r1 is not None:
			heho_records.append(r1)
		elif r2 is not None:
			other_records.append(r2)
	
	return (heho_records, other_records)

def create_family_vcf(v: tuple[VCFSmall, Family]) -> VCFFamily:
	return VCFFamily.create(v[0], v[1].samples())

def create_family_vcfs(orig_vcf: VCFSmall, large_families: list[KnownFamily],
										num_threads: int) -> list[VCFFamily]:
	v = [ (orig_vcf, family) for family in large_families ]
	if num_threads > 1:
		with Pool(num_threads) as pool:
			family_vcfs = pool.map(create_family_vcf, v)
	else:
		family_vcfs = list(map(create_family_vcf, v))
	return family_vcfs

def collect_same_parent_families(families: list[KnownFamily]
								) -> list[tuple[str, list[tuple[int, int]]]]:
	dic = defaultdict(list)		# { parent: [(family index, parent index)] }
	for i, family in enumerate(families):
		dic[family.mat].append((i, 0))
		dic[family.pat].append((i, 1))
	return sorted((p, v) for p, v in dic.items() if len(v) > 1)

def sort_records(rs: list[tuple[list[VCFHeteroHomoRecord],
								list[VCFImpFamilyRecord], int]]
							) -> list[list[tuple[VCFImpFamilyRecord, int]]]:
	indices1 = [ rs1[-1].index if rs1 else 0 for rs1, _, __ in rs ]
	indices2 = [ rs2[-1].index if rs2 else 0 for _, rs2, __ in rs ]
	max_index = max(indices1 + indices2)
	recordss: list[list[tuple[VCFImpFamilyRecord, int]]] = \
								[ [] for _ in range(max_index + 1) ]
	for rs1, rs2, parent_index in rs:
		for r1 in rs1:
			if r1.parents_wrong_type == 'Right':
				recordss[r1.index].append((r1, parent_index))
		for r2 in rs2:
			if r2.pair.is_homohomo():
				recordss[r2.index].append((r2, parent_index))
	return recordss

def modify_00x11_each(rs: list[tuple[list[VCFHeteroHomoRecord],
									 list[VCFImpFamilyRecord], int]]) -> None:
	if len(rs[0][0]) == 0 and len(rs[0][1]) == 0:
		return
	recordss = sort_records(rs)
	for records in recordss:
		if len(records) < 2:
			continue
		VCFImpFamilyRecord.modify_00x11(records)

# 0/0 x 1/1のrecordで別の家系の親になっているとき、
# 他のホモ×ホモやヘテロ×ホモとGenotypeが違うとき、修正する
def modify_00x11(heho_recordss: list[list[VCFHeteroHomoRecord]],
				 other_recordss: list[list[VCFImpFamilyRecord]],
				 families: list[KnownFamily]) -> None:
	# { parent: [(family index, parent index)] }
	fams = collect_same_parent_families(families)
	# やり方を変える
	# どの家系のどちらの親かだけ持っておけばよい
	for _, v in fams:
		rs = [ (heho_recordss[fam_index], other_recordss[fam_index], p_index)
												for fam_index, p_index in v ]
		modify_00x11_each(rs)

# divide into hetero x homo and other types
def divide_vcf_into_record_types(family_vcfs: list[VCFFamily],
									families: list[KnownFamily], option: Option
									) -> tuple[list[list[VCFHeteroHomoRecord]],
											   list[list[VCFImpFamilyRecord]]]:
	heho_recordss: list[list[VCFHeteroHomoRecord]] = []
	other_recordss: list[list[VCFImpFamilyRecord]] = []
	for vcf, family in zip(family_vcfs, families):
		heho_records, other_records = classify_records(vcf, family, option)
		heho_recordss.append(heho_records)
		other_recordss.append(other_records)
	return (heho_recordss, other_recordss)

def fill_vcf(dic_vcfs: dict[str, list[VCFHeteroHomo]],
						other_recordss: list[list[VCFImpFamilyRecord]],
						families: list[KnownFamily]) -> list[VCFFillable]:
	# 同じ家系でまとめる
	family_indices = { family.parents(): i
							for i, family in enumerate(families) }
	vcfss_heho: list[list[VCFHeteroHomo]] = [ [] for _ in families ]
	for vcfs in dic_vcfs.values():
		for vcf in vcfs:
			index = family_indices[vcf.parents()]
			vcfss_heho[index].append(vcf)
	
	# 家系ごとに残ったレコードをcorrectする
	return [ VCFFillable.fill(vcfs, other_records)
				for vcfs, other_records in zip(vcfss_heho, other_recordss) ]

def impute_hetero_homo(orig_vcf: VCFSmall, families: list[KnownFamily],
									geno_map: Map, op: Option
									) -> tuple[list[list[VCFHeteroHomo]],
											   list[list[VCFImpFamilyRecord]]]:
	# create a VCF for each large family
	family_vcfs = create_family_vcfs(orig_vcf, families, op.num_threads)
	# classify records
	heho_recordss, other_recordss = divide_vcf_into_record_types(family_vcfs,
																families, op)
	# 家系間で調整するから家系をまとめて処理しなければならない
	modify_00x11(heho_recordss, other_recordss, families)
	
	# 家系ごとに処理する
	vcfss: list[list[VCFHeteroHomo]] = []
	for records, others, family in zip(heho_recordss, other_recordss, families):
		vcfs, unused = VCFHeteroHomo.impute_vcfs(records, family.samples(),
															geno_map, orig_vcf)
		vcfss.append(vcfs)
		others.extend(unused)
	
	return (vcfss, other_recordss)

def impute(orig_vcf: VCFSmall, families: list[KnownFamily],
							  geno_map: Map, option: Option) -> VCFGeno:
	vcfss, other_recordss = impute_hetero_homo(orig_vcf, families,
															geno_map, option)
	
	# 同じ親でVCFをまとめる
	dic_vcfs: dict[str, list[VCFHeteroHomo]] = defaultdict(list)
	for vcfs, family in zip(vcfss, families):
		for vcf in vcfs:
			if len(vcf) == 0:
				# 家系で一つもVCFが無いと家系がわからなるので、
				# 空でも追加しておく
				dic_vcfs[family.mat].append(vcf)
				dic_vcfs[family.pat].append(vcf)
			elif vcf.is_mat_hetero():
				dic_vcfs[family.mat].append(vcf)
			else:
				dic_vcfs[family.pat].append(vcf)
	
	v = sorted(dic_vcfs.items())
	for parent, vcfs in v:
		VCFHeteroHomo.inverse_phases(vcfs)
	
	filled_vcfs = fill_vcf(dic_vcfs, other_recordss, families)
	print("%d large families have been imputed." % len(families))
	return VCFFillable.merge(filled_vcfs, orig_vcf.get_samples())

__all__ = ['impute_large_normal_families']
