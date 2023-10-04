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
from LargeFamily import correct_large_family_VCFs
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
		self.imputed_samples.update([ s for s in samples if s != '0' ])
	
	def clear(self):
		self.imputed_samples.clear()
	
	def is_imputed(self, sample: str) -> bool:
		return sample in self.imputed_samples
	
	def is_known(self, sample):
		return sample != '0'
	
	def is_unknown(self, sample):
		return sample == '0'
	
	def collet_references(self):
		return sorted(set(p for f in self.large_families for p in f.parents()))
	
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
							self.is_known(family.mat) and
							self.is_known(family.pat) and
							any(not self.is_imputed(prog)
											for prog in family.progenies) ]
		
		return [ Family(f.mat, f.pat, [ prog for prog in f.progenies
										if not self.is_imputed(prog) ])
															for f in families ]
	
	# 片親が補完されていて片親がunknownな家系
	def extract_phased_and_unknown_parents_family(self):
		families = [ family for family in self.small_families
					 if ((self.is_imputed(family.mat) and
					 						self.is_unknown(family.pat)) or
						  self.is_imputed(family.pat) and
						  					self.is_unknown(family.mat)) and
						  any(not self.is_imputed(prog)
											for prog in family.progenies) ]
		
		return [ Family(f.mat, f.pat, [ prog for prog in f.progenies
										if not self.is_imputed(prog) ])
															for f in families ]
	
	# 両親は補完されていないが子どもの一部が補完されている家系
	def extract_progenies_phased_families(self):
		families = [ family for family in self.small_families
						if (self.is_known(family.mat) or
									self.is_known(family.pat)) and
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
				samples.extend(s for s in family.samples() if self.is_known(s))
			elif self.is_unknown(family.mat) and self.is_unknown(family.pat):
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

def impute_vcf_by_parents(orig_vcf: VCFSmall, merged_vcf: VCFSmallBase,
								families: list[Family], gmap: Map) -> VCFSmall:
	vcfs: list[VCFSmallBase] = []
	for family in families:
		family_vcf = VCFHeteroHomoPP.impute_by_parents(
								orig_vcf, merged_vcf, family.samples(), gmap)
		vcf = family_vcf.remove_parents()
		vcfs.append(vcf)
	
	new_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	return new_vcf

def impute_vcf_by_parent(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
								families: list[Family], gmap: Map,
								sample_man: SampleManager,
								num_threads: int) -> VCFSmall:
	# collect not phased parents
	samples: list[str] = [ parent for family in families
								  for parent in family.parents()
								  if not sample_man.is_imputed(parent) ]
	
	# phase not phased parents
	parents_vcf = impute_iolated_samples(orig_vcf, merged_vcf, sample_man,
													samples, gmap, num_threads)
	
	# merge vcfs
	new_merged_vcf = VCFSmall.join([merged_vcf, parents_vcf], orig_vcf.samples)
	
	# impute progenies
	new_vcf = impute_vcf_by_parents(orig_vcf, new_merged_vcf, families, gmap)
	
	# join
	return VCFSmall.join([parents_vcf, new_vcf], orig_vcf.samples)

def impute_one_parent_vcf(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
								families: list[Family], gmap: Map,
								sample_man: SampleManager,
								num_threads: int) -> VCFSmall:
	references = sample_man.collet_references()
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

def impute_vcf_by_progenies(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
									families: list[Family], gmap: Map,
									sample_man: SampleManager) -> VCFSmall:
	references = sample_man.collet_references()
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
			vcf2 = impute_vcf_by_parents(orig_vcf, pp_vcf, families, gmap)
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

def impute_iolated_samples(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
									sample_man: SampleManager,
									samples: list[str],
									gmap: Map, num_threads: int) -> VCFSmall:
	references = sample_man.collet_references()
	# あとでマルチプロセス化するためにphasingすべきsample分割する
	vcfs = VCFIsolated.create(orig_vcf, merged_vcf,
								samples, references, gmap, num_threads)
	new_vcfs: list[VCFSmallBase] = []
	for vcf in vcfs:
		vcf.impute()
		new_vcfs.append(vcf.extract_isolated_samples())
	
	new_vcf = VCFSmall.join(new_vcfs, orig_vcf.samples)
	return new_vcf

def impute_vcf_chr(orig_vcf: VCFSmall, sample_man: SampleManager,
						geno_map: Map, option: Option) -> VCFSmall:
	print('chr: %s %d records' % (orig_vcf.records[0].chrom(), len(orig_vcf)))
	sys.stdout.flush()
	merged_vcf = correct_large_family_VCFs(orig_vcf, sample_man.large_families,
														geno_map, option)
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
												   families, geno_map,
												   sample_man,
												   option.num_threads)
			merged_vcf = VCFSmall.join([merged_vcf, new_imputed_vcf],
														orig_vcf.samples)
			sample_man.set(new_imputed_vcf.samples)
			continue
		
		families = sample_man.extract_phased_and_unknown_parents_family()
		if families:
			new_imputed_vcf = impute_one_parent_vcf(orig_vcf, merged_vcf,
													families, geno_map,
													sample_man,
													option.num_threads)
			merged_vcf = VCFSmall.join([merged_vcf, new_imputed_vcf],
														orig_vcf.samples)
			sample_man.set(new_imputed_vcf.samples)
			continue
		
		families = sample_man.extract_progenies_phased_families()
		if families:
			new_imputed_vcf = impute_vcf_by_progenies(orig_vcf, merged_vcf, 
														families, geno_map,
														sample_man)
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
		vcfs: list[VCFSmallBase] = [merged_vcf, new_imputed_vcf]
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

def print_info(option: Option):
	print("input VCF : %s" % option.path_VCF, file=sys.stderr)
	print("pedigree : %s" % option.path_ped, file=sys.stderr)
	print("output VCF : %s" % option.path_out, file=sys.stderr)

def impute_vcf(option: Option):
	print_info(option)
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
