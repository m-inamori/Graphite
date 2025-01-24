from __future__ import annotations

# coding: utf-8
# SmallFamily.py

from itertools import *
from collections import defaultdict, Counter
import random
import csv
import sys
from multiprocessing import Pool

from typing import Dict, List, Tuple, Set, Optional

from VCF import *
from VCFFamily import VCFFamily, VCFFamilyBase, VCFFamilyRecord
from pedigree import PedigreeTable, Family
from VCFHeteroHomoPP import *
import OnePhasedFamily
from VCFOneParentPhased import VCFOneParentPhased
import NoPhasedFamily 
from VCFProgenyPhased import VCFProgenyPhased
from VCFIsolated import VCFIsolated
from SampleManager import *
from Map import *
from option import *
from common import *


def impute_vcf_by_parents_core(orig_vcf: VCFSmall, merged_vcf: VCFSmallBase,
							families: list[KnownFamily], gmap: Map) -> VCFSmall:
	vcfs: list[VCFSmallBase] = []
	for family in families:
		family_vcf = VCFHeteroHomoPP.impute_by_parents(
								orig_vcf, merged_vcf, family.samples(), gmap)
		vcf = family_vcf.remove_parents()
		vcfs.append(vcf)
	
	new_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	return new_vcf

def impute_vcf_by_parents(orig_vcf: VCFSmall, merged_vcf: VCFSmallBase,
							sample_man: SampleManager,
							gmap: Map, num_threads: int) -> Optional[VCFSmall]:
	families = sample_man.extract_small_families()
	if not families:
		return None
	
	new_imputed_vcf = impute_vcf_by_parents_core(orig_vcf, merged_vcf,
															families, gmap)
	merged_vcf = VCFSmall.join([merged_vcf, new_imputed_vcf], orig_vcf.samples)
	sample_man.set(new_imputed_vcf.samples)
	return merged_vcf

def impute_vcf_by_parent(orig_vcf: VCFSmall, imputed_vcf: VCFSmall,
										ref_haps: list[list[int]], gmap: Map,
										sample_man: SampleManager,
										num_threads: int) -> Optional[VCFSmall]:
	families = sample_man.extract_single_parent_phased_families()
	if not families:
		return None
	
	# families have been already selected
	# in which one parent has been imputed and one parent has not been imputed
	# collect not phased parents
	samples: list[str] = [ parent for family in families
								  for parent in family.parents()
								  if not sample_man.is_imputed(parent) ]
	
	vcf = OnePhasedFamily.impute_by_parent(orig_vcf, imputed_vcf,
											ref_haps, families, samples, gmap)
	merged_vcf = VCFSmall.join([imputed_vcf, vcf], orig_vcf.samples)
	sample_man.set(vcf.get_samples())
	return merged_vcf

def impute_vcf_by_non_imputed_parents(orig_vcf: VCFSmall, imputed_vcf: VCFSmall,
							ref_haps: list[list[int]], gmap: Map,
							sample_man: SampleManager) -> Optional[VCFSmall]:
	families = sample_man.extract_no_parent_phased_families()
	if not families:
		return None
	
	vcf = NoPhasedFamily.impute(orig_vcf, imputed_vcf, ref_haps, families, gmap)
	merged_vcf = VCFSmall.join([imputed_vcf, vcf], orig_vcf.samples)
	sample_man.set(vcf.get_samples())
	return merged_vcf

def impute_one_parent_vcf_core(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
										families: list[Family], gmap: Map,
										sample_man: SampleManager,
										num_threads: int) -> VCFSmall:
	references = sample_man.collect_reference()
	ref_vcf = merged_vcf.extract_samples(references)
	
	vcfs: list[VCFSmallBase] = []
	for family in families:
		is_mat_phased = sample_man.is_imputed(family.mat)
		opp_vcf = VCFOneParentPhased.create(family.samples(), is_mat_phased,
											merged_vcf, orig_vcf, gmap, ref_vcf)
		opp_vcf.impute()
		vcfs.append(opp_vcf)
	
	samples = []	# 今回phasingしたsampleを集める
	for family in families:
		ss = [ s for s in family.samples()
					if not sample_man.is_imputed(s) and s != '0' ]
		samples.extend(ss)
		sample_man.set(ss)
	
	return VCFSmall.join(vcfs, samples)

def impute_one_parent_vcf(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
										gmap: Map, sample_man: SampleManager,
										num_threads: int) -> Optional[VCFSmall]:
	families = sample_man.extract_phased_and_unknown_parents_family()
	if not families:
		return None
	
	new_imputed_vcf = impute_one_parent_vcf_core(orig_vcf, merged_vcf, families,
												 gmap, sample_man, num_threads)
	merged_vcf = VCFSmall.join([merged_vcf, new_imputed_vcf],
												orig_vcf.samples)
	sample_man.set(new_imputed_vcf.samples)
	return merged_vcf

def impute_vcf_by_progenies_core(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
											families: list[Family], gmap: Map,
											sample_man: SampleManager,
											num_threads: int) -> VCFSmall:
	references = sample_man.collect_reference()
	ref_vcf = merged_vcf.extract_samples(references)
	vcfs: list[VCFSmallBase] = []
	for family in families:
		ppi = next(i for i, sample in enumerate(family.samples())
								if sample_man.is_imputed(sample))
		pp_vcf = VCFProgenyPhased.create(orig_vcf, merged_vcf,
									family.samples(), ppi, gmap, ref_vcf)
		pp_vcf.impute()		# impute only parents
		if pp_vcf.num_progenies() == 1:
			vcf4: VCFSmallBase = pp_vcf	# all imputed
		elif family.mat != '0' and family.pat != '0':
			vcf2 = impute_vcf_by_parents(orig_vcf, pp_vcf,
											sample_man, gmap, num_threads)
			if vcf2 is None:
				return merged_vcf	# not comes here
			vcf4 = VCFSmall.join([pp_vcf, vcf2], family.samples())
		else:
			is_mat_phased = family.mat != '0'
			vcf3 = VCFOneParentPhased.create(family.samples(), is_mat_phased,
												pp_vcf, orig_vcf, gmap, ref_vcf)
			vcf3.impute()
			vcf4 = VCFSmall.join([pp_vcf, vcf3], family.samples())
		vcfs.append(vcf4)
	
	samples = []	# collect phased samples
	for family in families:
		ss = [ s for s in family.samples()
							if not sample_man.is_imputed(s) ]
		samples.extend(ss)
		sample_man.set(ss)
	
	return VCFSmall.join(vcfs, samples)

def impute_vcf_by_progenies(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
									gmap: Map, sample_man: SampleManager,
									num_threads: int) -> Optional[VCFSmall]:
	families = sample_man.extract_progenies_phased_families()
	if not families:
		return None
	
	new_imputed_vcf = impute_vcf_by_progenies_core(orig_vcf, merged_vcf,
													families, gmap,
													sample_man, num_threads)
	merged_vcf = VCFSmall.join([merged_vcf, new_imputed_vcf],
												orig_vcf.samples)
	sample_man.set(new_imputed_vcf.samples)
	return merged_vcf

def extract_haplotypes(phased_vcf: VCFSmall,
						sample_man: SampleManager) -> list[list[int]]:
	reference: list[str] = sample_man.collect_reference()
	ref_columns = phased_vcf.extract_columns(reference)
	records: list[VCFRecord] = []
	for record in phased_vcf.records:
		v = record.v[:9] + [ record.v[c] for c in ref_columns ]
		new_record = VCFRecord(v, reference)
		records.append(new_record)
	
	N = len(reference)
	M = len(records)
	
	def geno(h: int, i: int) -> int:
		return int(records[i].v[(h>>1)+9][(h&1)<<1])
	
	gts = [ [ geno(h, i) for i in range(M) ] for h in range(N*2) ]
	
	K = min(10, M)
	# ローリングハッシュ
	a: list[list[int]] = [ [0] * (M-K+1) for _ in range(N*2) ]
	for h in range(N*2):
		n = 0
		for i in range(M):
			gt = gts[h][i]
			if i < K - 1:
				n = (n << 1) | gt
			else:
				n = ((n << 1) & ((1 << 10) - 1)) | gt
				a[h][i-K+1] = n
	
	# 同じハッシュ値があったらどのハプロタイプと同じか調べる
	b: list[list[int]] = [ [0] * (M-K+1) for _ in range(N*2) ]
	for i in range(M-K+1):
		for h in range(N*2):
			for l in range(h):
				if a[l][i] == a[h][i]:
					b[h][i] = l
					break
			else:
				b[h][i] = h
	
	# Genotypeに戻ってできるだけ自分以外になるようにする
	for k in range(1, N*2):
		i = 0
		c: list[tuple[int, int, int]] = []	# [(index, first, last)]
		for h, v1 in groupby(b[k]):
			next_i = i + sum(1 for _ in v1)
			c.append((h, i, next_i))
			i = next_i
		
		# 自分以外の前後を可能なら拡張する
		for j in range(len(c)):
			h0, first0, last0 = c[j]
			if h0 != k:
				continue
			if j != 0:
				# 前を見る
				h1, first1, last1 = c[j-1]
				for i in range(first0, last0):
					if gts[h1][i] != gts[k][i]:
						break
					b[k][i] = h1
			if j != len(c) - 1:
				# 後ろを見る
				h2, first2, last2 = c[j+1]
				for i in range(last0-1, first0, -1):
					if gts[h2][i] != gts[k][i]:
						break
					b[k][i] = h2
	
	# 2つのハプロタイプが長さの1割以下しかないサンプルは捨てる
	counter = Counter(g for v in b for g in v)
	indices = [ i for i in range(N*2) if counter[i] * 10 >= M ]
	ref_gts = [ gts[i] for i in indices ]
	return ref_gts

def impute_small_family_VCFs(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
								geno_map: Map, sample_man: SampleManager,
								num_threads: int) -> VCFSmall:
	ref_haps = extract_haplotypes(merged_vcf, sample_man)
	while True:
		# 両親が補完されているが子どもが少ない家系を補完する
		new_merged_vcf1 = impute_vcf_by_parents(orig_vcf, merged_vcf,
											sample_man, geno_map, num_threads)
		if new_merged_vcf1 is not None:
			merged_vcf = new_merged_vcf1
			continue
		
		# 片親が補完されている家系を補完する
		new_merged_vcf2 = impute_vcf_by_parent(orig_vcf, merged_vcf,
												ref_haps, geno_map,
												sample_man, num_threads)
		if new_merged_vcf2 is not None:
			merged_vcf = new_merged_vcf2
			continue
		
		# 両親が補完されていない家系を補完する
		new_merged_vcf3 = impute_vcf_by_non_imputed_parents(orig_vcf,
														merged_vcf, ref_haps,
														geno_map, sample_man)
		if new_merged_vcf3 is not None:
			merged_vcf = new_merged_vcf3
			continue
		
		new_merged_vcf4 = impute_one_parent_vcf(orig_vcf, merged_vcf,
											geno_map, sample_man, num_threads)
		if new_merged_vcf4 is not None:
			merged_vcf = new_merged_vcf4
			continue
		
		# Impute families whose progenies have been imputed
		new_merged_vcf5 = impute_vcf_by_progenies(orig_vcf, merged_vcf,
											geno_map, sample_man, num_threads);
		if new_merged_vcf5 is not None:
			merged_vcf = new_merged_vcf5
			continue
		
		break
	
	return merged_vcf

def impute_isolated_samples(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
									sample_man: SampleManager,
									samples: list[str],
									gmap: Map, num_threads: int) -> VCFSmall:
	reference = sample_man.collect_reference()
	# あとでマルチプロセス化するためにphasingすべきsample分割する
	vcfs = VCFIsolated.create(orig_vcf, merged_vcf,
								samples, reference, gmap, num_threads)
	new_vcfs: list[VCFSmallBase] = []
	for vcf in vcfs:
		vcf.impute()
		new_vcfs.append(vcf.extract_isolated_samples())
	
	new_vcf = VCFSmall.join(new_vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute_small_family_VCFs', 'impute_isolated_samples']
