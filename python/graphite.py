from __future__ import annotations

# coding: utf-8
# graphite.py
# フルでphasingする

from itertools import *
from collections import defaultdict, Counter
import random
import csv
import sys
from multiprocessing import Pool

from typing import Dict, List, Tuple, Set, IO

from VCF import *
from VCFFamily import VCFFamily, VCFFamilyBase
from pedigree import PedigreeTable, Family
from VCFHeteroHomo import *
from VCFHomoHomo import VCFHomoHomoRecord
from VCFHeteroHeteroLite import VCFHeteroHeteroLiteRecord
from VCFJunkRecord import *
from VCFFillable import *
from VCFHeteroHomoPP import *
from VCFOneParentPhased import *
from VCFProgenyPhased import *
from VCFIsolated import *
from Map import *
import ClassifyRecord as CR
from option import *
from common import *


#################### library ####################


#################### SampleManager ####################

class SampleManager:
	def __init__(self, ped: PedigreeTable, large_families: list[Family],
									small_families: list[Family], lower_p: int):
		self.ped: PedigreeTable				= ped
		self.large_families: list[Family]	= large_families
		self.small_families: list[Family]	= small_families
		self.lower_progs: int				= lower_p
		self.imputed_samples: set[str]		= set()
	
	def set(self, samples: list[str]):
		self.imputed_samples.update(samples)
	
	def clear(self):
		self.imputed_samples.clear()
	
	def is_imputed(self, sample: str) -> bool:
		return sample in self.imputed_samples
	
	# 補完されていないが両親は補完されている家系
	def extract_small_families(self):
		families = [ family for family in self.small_families
						if self.is_imputed(family.mat) and
							self.is_imputed(family.pat) and
							any(not self.is_imputed(prog)
											for prog in family.progenies) ]
		
		return [ Family(f.mat, f.pat, [ prog for prog in f.progenies
										if not self.is_imputed(prog) ])
															for f in families ]
	
	# 補完されていないが片方の親だけ補完されている家系
	def extract_single_parent_phased_families(self):
		families = [ family for family in self.small_families
						if (self.is_imputed(family.mat) or
							self.is_imputed(family.pat)) and
							any(not self.is_imputed(prog)
											for prog in family.progenies) ]
		
		return [ Family(f.mat, f.pat, [ prog for prog in f.progenies
										if not self.is_imputed(prog) ])
															for f in families ]
	
	# 両親は補完されていないが子どもの一部が補完されている家系
	def extract_progenies_phased_families(self):
		families = [ family for family in self.small_families
						if (family.mat != '0' or family.pat != '0') and
							not self.is_imputed(family.mat) and
							not self.is_imputed(family.pat) and
							any(self.is_imputed(prog)
										for prog in family.progenies) ]
		return families
	
	def extract_isolated_samples(self):
		# 繋がっているサンプルがあっても、
		# 家系の全サンプルがphasingされていないなら孤立とみなす
		samples = []
		for family in self.small_families:
			if all(not self.is_imputed(s) for s in family.samples()):
				samples.extend(s for s in family.samples() if s != '0')
			elif family.mat == '0' and family.pat == '0':
				# 親が両方とも不明なら、phasingされていないサンプルはOK
				samples.extend(s for s in family.progenies
											if not self.is_imputed(s))
		return samples
	
	def display_info(self, out: IO):
		print("%d samples" % len(self.ped), file=out)
		
		if len(self.large_families) == 1:
			print("1 large family", end='', file=out)
		else:
			print("%d large families" % len(self.large_families), end='',
														file=out)
		print(" (number of progenies >= %d)" % self.lower_progs, file=out)
		
		if len(self.small_families) == 1:
			print("1 small family", file=out)
		else:
			print("%d small families" % len(self.small_families), file=out)
	
	@staticmethod
	def make_families(ped: PedigreeTable,
						set_samples: Set[str]) -> list[Family]:
		# 不明扱いを後ですることにより、両親が同じサンプルをひとまとめにする
		dic = classify((prog.parents(), prog.name) for prog in ped.table)
		families = []
		for (mat, pat), progs in dic.items():
			filtered_progs = [ prog for prog in progs if prog in set_samples ]
			if not filtered_progs:
				continue
			# VCFに無い親は不明扱い
			mat_mod = mat if mat in set_samples else '0'
			pat_mod = pat if pat in set_samples else '0'
			family = Family(mat_mod, pat_mod, filtered_progs)
			families.append(family)
		families.sort(key=lambda f: f.parents())
		return families
	
	@staticmethod
	def create(path_ped: str, samples: list[str],
				lower_progs: int, family_indices: list[int]) -> SampleManager:
		ped_ = PedigreeTable.read(path_ped)
		ped = ped_.limit_samples(samples)
		
		set_samples = set(samples)
		families = SampleManager.make_families(ped, set_samples)
		if family_indices:	# debug用にFamilyを絞って問題を小さくする
			families = [ families[i] for i in family_indices ]
		
		large_families = []
		small_families = []
		for f in families:
			if (f.num_progenies() >= lower_progs and
							f.mat in set_samples and f.pat in set_samples):
				large_families.append(f)
			else:
				small_families.append(f)
		return SampleManager(ped, large_families, small_families, lower_progs)


#################### process ####################

# FamilyごとにVCFHeteroHomoを作って親ごとに格納する
def make_VCFHeteroHomo(records: list[VCFHeteroHomoRecord],
						family: Family, header: list[list[str]], geno_map: Map
								) -> tuple[VCFHeteroHomo, VCFHeteroHomo,
												list[VCFHeteroHomoRecord]]:
	unused_records: list[VCFHeteroHomoRecord] = []
	samples = family.samples()
	new_header = header[:-1] + [header[-1][:9] + samples]
	heho_mat_records = [ r for r in records
							if r.is_imputable() and r.is_mat_hetero() ]
	heho_pat_records = [ r for r in records
							if r.is_imputable() and not r.is_mat_hetero() ]
	unused_records = [ r for r in records if not r.is_imputable() ]
	vcf_mat = VCFHeteroHomo(new_header, heho_mat_records, geno_map)
	vcf_pat = VCFHeteroHomo(new_header, heho_pat_records, geno_map)
	return (vcf_mat, vcf_pat, unused_records)

def impute_each_parent(v: tuple[list[VCFHeteroHomo], Option]
					) -> tuple[list[VCFHeteroHomo], list[VCFHeteroHomoRecord]]:
	vcfs_heho, option = v
	return VCFHeteroHomo.impute_vcfs(vcfs_heho, option)

def impute_hetero_homo_parellel(vcfs: Dict[str,list[VCFHeteroHomo]],
							option: Option) -> list[tuple[list[VCFHeteroHomo],
													list[VCFHeteroHomoRecord]]]:
	v = [ (vcfs_heho, option) for parent, vcfs_heho in vcfs.items() ]
	num_threads = min(option.num_threads, len(v))
	if num_threads > 1 and len(v[0][0][0]) >= 10:
		with Pool(num_threads) as pool:
			return pool.map(impute_each_parent, v)
	else:
		return list(map(impute_each_parent, v))

def impute_hetero_homo_core(records: Dict[Family, list[VCFHeteroHomoRecord]],
						header: list[list[str]], geno_map: Map, option: Option
							) -> tuple[Dict[Family, list[VCFHeteroHomo]],
									   Dict[Family, list[VCFHeteroHomoRecord]]]:
	# あとでFamilyを回復するために必要
	parents_to_family = { family.parents(): family
										for family in records.keys() }
	vcfs: Dict[str,list[VCFHeteroHomo]] = defaultdict(list)
	unused_records: Dict[Family, list[VCFHeteroHomoRecord]] = defaultdict(list)
	for family, heho_records in records.items():
		vcf_mat, vcf_pat, unused_records1 = make_VCFHeteroHomo(heho_records,
													family, header, geno_map)
		vcfs[family.mat].append(vcf_mat)
		vcfs[family.pat].append(vcf_pat)
		unused_records[family].extend(unused_records1)
	
	w = impute_hetero_homo_parellel(vcfs, option)
	
	imputed_vcfs: Dict[Family, list[VCFHeteroHomo]] = defaultdict(list)
	for pair in w:
		imputed_vcfs_heho, unused_records1 = pair
		for vcf in imputed_vcfs_heho:
			family = parents_to_family[vcf.parents()]
			imputed_vcfs[family].append(vcf)
		
		for record in unused_records1:
			record.enable_modification()
			parents = (record.samples[0], record.samples[1])
			family = parents_to_family[parents]
			unused_records[family].append(record)
	return (imputed_vcfs, unused_records)

Item = Tuple[List[VCFHeteroHomo], List[VCFImpFamilyRecord]]

def fill_vcf(v: Item) -> VCFFillable:
	vcfs, records = v
	return VCFFillable.fill(vcfs, records)

# HeteroHomoだけ別にする
# このあとHeteroHomoだけ補完するから
# その他はVCFFillableにした後補完する
def classify_records(vcf: VCFSmall, families: list[Family], p: float
							) -> tuple[Dict[Family, list[VCFHeteroHomoRecord]],
									   Dict[Family, list[VCFImpFamilyRecord]]]:
	heho_records: Dict[Family, list[VCFHeteroHomoRecord]] = \
									{ family: [] for family in families }
	other_records: Dict[Family, list[VCFImpFamilyRecord]] = \
									{ family: [] for family in families }
	for family in families:
		vcf_family = VCFFamily.create(vcf, family.samples())
		samples = family.samples()
		td = CR.get_typedeterminer(family.num_progenies(), p)
		for i, record in enumerate(vcf_family.records):
			wrong_type, pair = CR.classify_record(record, td)
			v = record.v
			if pair == -1:	# 候補が無い
				record = VCFJunkRecord(v, samples, i, wrong_type)
				other_records[family].append(record)
			elif pair == 0 or pair == 3 or pair == 5:
				record_ = VCFHomoHomoRecord(v, samples, i, wrong_type, pair)
				record = record_.impute()
				other_records[family].append(record)
			elif pair == 1 or pair == 4:		# 0/0 x 0/1 or 0/1 x 1/1
				record = VCFHeteroHomoRecord(v, samples, i, wrong_type, pair)
				heho_records[family].append(record)
			else:		# 0/1 x 0/1
				record = VCFHeteroHeteroLiteRecord(v, samples, i, wrong_type)
				other_records[family].append(record)
	
	return (heho_records, other_records)

def sort_records(heho_records: Dict[Family, list[VCFHeteroHomoRecord]],
				 other_records: Dict[Family, list[VCFImpFamilyRecord]]
		 									) -> list[list[VCFImpFamilyRecord]]:
	is1 = [ r.index for rs in heho_records.values() for r in rs ]
	is2 = [ r.index for rs in other_records.values() for r in rs ]
	max_index = max(is1 + is2)
	records: list[list[VCFImpFamilyRecord]] = \
							[ [] for _ in range(max_index + 1) ]
	for rs1 in heho_records.values():
		for r1 in rs1:
			if r1.parents_wrong_type == 'Right':
				records[r1.index].append(r1)
	for rs2 in other_records.values():
		for r2 in rs2:
			if r2.pair in (0, 3, 5):		# Homo x Homo
				records[r2.index].append(r2)
	return records

# 0/0 x 1/1のrecordで別の家系の親になっているとき、
# 他のホモ×ホモやヘテロ×ホモとGenotypeが違うとき、修正する
def modify_00x11(heho_records: Dict[Family, list[VCFHeteroHomoRecord]],
				 other_records: Dict[Family, list[VCFImpFamilyRecord]]):
	records = sort_records(heho_records, other_records)
	
	for rs in records:
		if len(rs) < 2:
			continue
		VCFImpFamilyRecord.modify_00x11(rs)

def impute_hetero_homo(orig_vcf: VCFSmall, sample_man: SampleManager,
										geno_map: Map, option: Option
							) -> tuple[Dict[Family, list[VCFHeteroHomo]],
									   Dict[Family, list[VCFImpFamilyRecord]]]:
	for f in sample_man.large_families:
		CR.prepare(len(f.samples()) - 2, option.ratio)
	heho_records, other_records = classify_records(orig_vcf,
										sample_man.large_families, option.ratio)
	modify_00x11(heho_records, other_records)
	
	# HeteroHomoとそれ以外を分けて、HeteroHomoを補完する
	# 使われなかったHeteroHomoはvillに回す
	imputed_vcfs, unused_records = impute_hetero_homo_core(heho_records,
											orig_vcf.header, geno_map, option)
	
	for family, records in unused_records.items():
		other_records[family].extend(records)
		other_records[family].sort(key=lambda r: r.index)
	
	return (imputed_vcfs, other_records)

def fill_and_merge_vcf(imputed_vcfs: Dict[Family, list[VCFHeteroHomo]],
						other_records: Dict[Family, list[VCFImpFamilyRecord]],
						samples: list[str], option: Option) -> VCFSmall:
	def fill(items: list[Item]) -> list[VCFFillable]:
		num_threads = min(len(items), option.num_threads)
		vcfs, records = items[0]
		num_records = sum(map(len, vcfs)) + len(records)
		if num_threads > 1 and num_records >= 100:
			with Pool(option.num_threads) as pool:
				return pool.map(fill_vcf, items)
		else:
			return list(map(fill_vcf, items))
	
	items = [ (vcfs, other_records[family])
							for family, vcfs in imputed_vcfs.items() ]
	# C++と挙動を合わせるため
	items.sort(key=lambda v: tuple(v[0][0].samples[:2]))
	vcfs_complemented: list[VCFFillable] = fill(items)
	return VCFFillable.merge(vcfs_complemented, samples)

def impute_vcf_by_parents(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
								families: list[Family], gmap: Map) -> VCFSmall:
	vcfs: list[VCFSmall] = []
	for family in families:
		family_vcf = VCFHeteroHomoPP.impute_by_parents(
								orig_vcf, merged_vcf, family.samples(), gmap)
		vcf = family_vcf.remove_parents()
		vcfs.append(vcf)
	
	new_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	return new_vcf

def impute_vcf_by_parent(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
								families: list[Family], gmap: Map,
								sample_man: SampleManager) -> VCFSmall:
	vcfs: list[VCFFamily] = []
	for family in families:
		is_mat_phased = sample_man.is_imputed(family.mat)
		opp_vcf = VCFOneParentPhased.impute_by_parent(
										orig_vcf, merged_vcf,
										family.samples(), is_mat_phased, gmap)
		family_vcf = VCFFamily.convert(opp_vcf)
		vcfs.append(family_vcf)
	
	samples = []	# 今回phasingした親と子どもを集める
	for family in families:
		samples.extend(family.progenies)
		if family.mat == '0' or family.pat == '0':
			continue
		elif sample_man.is_imputed(family.mat):
			samples.append(family.pat)
		else:
			samples.append(family.mat)
	
	merged_records: list[VCFRecord] = []
	for i in range(len(vcfs[0])):
		merged_records.append(VCFOneParentPhased.merge_records(
													vcfs, i, samples))
	new_header = orig_vcf.trim_header(samples)
	new_vcf = VCFSmall(new_header, merged_records)
	return new_vcf

def impute_vcf_by_progenies(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
									families: list[Family],
									sample_man: SampleManager) -> VCFSmall:
	vcfs: list[VCFSmall] = []
	for family in families:
		ppi = [ i for i, sample in enumerate(family.samples())
								if sample_man.is_imputed(sample) ]
		pp_vcf = VCFProgenyPhased.impute_by_progeny(
								orig_vcf, merged_vcf, family.samples(), ppi)
		family_vcf = VCFSmall.convert(pp_vcf)
		vcfs.append(family_vcf)
	
	new_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	return new_vcf

def impute_iolated_samples(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
									sample_man: SampleManager,
									samples: list[str],
									gmap: Map, num_threads: int) -> VCFSmall:
	references = sorted(set(p for f in sample_man.large_families
												for p in f.parents()))
	# あとでマルチプロセス化するためにphasingすべきsample分割する
	vcfs = VCFIsolated.create(orig_vcf, merged_vcf,
								samples, references, gmap, num_threads)
	new_vcfs: list[VCFSmall] = []
	for vcf in vcfs:
		vcf.impute()
		new_vcfs.append(vcf.extract_isolated_samples())
	
	new_vcf = VCFSmall.join(new_vcfs, orig_vcf.samples)
	return new_vcf

def impute_vcf_chr(orig_vcf: VCFSmall, sample_man: SampleManager,
						geno_map: Map, option: Option) -> VCFSmall:
	print("chr : ", orig_vcf.records[0].chrom())
	imputed_vcfs, other_records = impute_hetero_homo(orig_vcf, sample_man,
															geno_map, option)
	merged_vcf = fill_and_merge_vcf(imputed_vcfs,
									other_records, orig_vcf.samples, option)
	if option.only_large_families:
		return merged_vcf
	
	# 補完できる家系がなくなるまで繰り返す
	sample_man.set(merged_vcf.samples)
	while True:
		# 両親が補完されているが子どもが少ない家系を補完する
		families = sample_man.extract_small_families()
		if families:
			new_imputed_vcf = impute_vcf_by_parents(orig_vcf, merged_vcf,
														families, geno_map)
			merged_vcf = VCFSmall.join([merged_vcf, new_imputed_vcf],
														orig_vcf.samples)
			sample_man.set(new_imputed_vcf.samples)
			continue
		
		# 片親が補完されている家系を補完する
		families = sample_man.extract_single_parent_phased_families()
		if families:
			new_imputed_vcf = impute_vcf_by_parent(orig_vcf, merged_vcf, 
												families, geno_map, sample_man)
			merged_vcf = VCFSmall.join([merged_vcf, new_imputed_vcf],
														orig_vcf.samples)
			sample_man.set(new_imputed_vcf.samples)
			continue
		
		families = sample_man.extract_progenies_phased_families()
		if families:
			new_imputed_vcf = impute_vcf_by_progenies(orig_vcf, merged_vcf, 
														families, sample_man)
			merged_vcf = VCFSmall.join([merged_vcf, new_imputed_vcf],
														orig_vcf.samples)
			sample_man.set(new_imputed_vcf.samples)
			continue
		else:
			break
	
	# 最後に孤立したサンプルを補完する
	samples = sample_man.extract_isolated_samples()
	if samples:
		new_imputed_vcf = impute_iolated_samples(orig_vcf, merged_vcf,
							sample_man, samples, geno_map, option.num_threads)
		vcfs = [merged_vcf, new_imputed_vcf]
		merged_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	
	sample_man.clear()
	return merged_vcf

def chroms_efficients(option: Option) -> Iterator[bool]:
	if not option.chroms:
		return repeat(True)
	
	upper = max(option.chroms)
	v: list[bool] = [False] * (upper + 1)
	for i in option.chroms:
		v[i] = True
	return (b for b in v)

def impute_vcf(option: Option):
	vcf = VCFHuge.read(option.path_VCF)
	
	geno_map = Map.read(option.path_map)
	geno_map.display_info(sys.stderr)
	
	sample_man = SampleManager.create(option.path_ped, vcf.samples,
										option.lower_progs, option.families)
	sample_man.display_info(sys.stderr)
	
	iter = chroms_efficients(option)
	first = True
	for b, vcf_chr, gmap in zip(iter, vcf.divide_into_chromosomes(),
													geno_map.iter_chr_maps()):
		if not b:
			continue
		vcf_imputed = impute_vcf_chr(vcf_chr, sample_man, gmap, option)
		with open(option.path_out, 'w' if first else 'a') as out:
			vcf_imputed.write(out, with_header=first)
		first = False


#################### main ####################

option = Option.create(sys.argv)
if option is None:
	Option.usage()
	exit(1)

impute_vcf(option)
