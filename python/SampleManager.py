from __future__ import annotations

# coding: utf-8
# SampleManager.py

from typing import Dict, List, Tuple, Set, IO

from pedigree import PedigreeTable, Family
from common import classify


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
