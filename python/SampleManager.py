from __future__ import annotations

# coding: utf-8
# SampleManager.py

from typing import Set, TextIO, Optional

from pedigree import PedigreeTable, Family
from KnownFamily import KnownFamily
from common import classify
from exception_with_code import *
import error_codes


#################### SampleException ####################

class SampleException(ExceptionWithCode):
	def __init__(self) -> None:
		super().__init__('error : large family not found')
	
	def get_error_code(self) -> error_codes.Type:
		return error_codes.Type.LARGE_FAMILY_NOT_FOUND


#################### SampleManager ####################

class SampleManager:
	def __init__(self, ped: PedigreeTable, large_families: list[KnownFamily],
							large_self_families: list[KnownFamily],
							small_families: list[KnownFamily], lower_p: int):
		self.ped: PedigreeTable				= ped
		self.large_families: list[KnownFamily]	= large_families
		self.large_self_families: list[KnownFamily]	= large_self_families
		self.small_families: list[KnownFamily]	= small_families
		self.lower_progs: int				= lower_p
		self.imputed_samples: set[str]		= set()
	
	def set(self, samples: list[str]) -> None:
		self.imputed_samples.update([ s for s in samples if s != '0' ])
	
	def clear(self) -> None:
		self.imputed_samples.clear()
	
	def is_imputed(self, sample: str) -> bool:
		return sample in self.imputed_samples
	
	def is_known(self, sample: str) -> bool:
		return sample != '0'
	
	def collect_reference(self) -> list[str]:
		s = set(p for f in self.large_families
							for p in f.known_parents())
		self_families = self.extract_self_parent_imputed_families()
		for family in self_families:
			s.add(family.mat)
		return sorted(s)
	
	def extract_unimputed_progenies(self, f: Family) -> list[str]:
		return [ prog for prog in f.progenies if not self.is_imputed(prog) ]
	
	# 補完されていないが両親は補完されている家系
	def extract_both_imputed_families(self) -> list[KnownFamily]:
		families = [ family for family in self.small_families
						if not family.is_self() and
							self.is_imputed(family.mat) and
							self.is_imputed(family.pat) and
							any(not self.is_imputed(prog)
											for prog in family.progenies) ]
		
		return [ KnownFamily(f.mat, f.pat, f.mat_known, f.pat_known,
									self.extract_unimputed_progenies(f))
															for f in families ]
	
	# 補完されていないが片方の親だけ補完されている家系
	def extract_imputed_and_known_families(self) -> list[Family]:
		families = [ family for family in self.small_families
						if not family.is_self() and
							(self.is_imputed(family.mat) ^
							self.is_imputed(family.pat)) and
							family.mat_known and family.pat_known and
							any(not self.is_imputed(prog)
											for prog in family.progenies) ]
		
		return [ Family(f.mat, f.pat, self.extract_unimputed_progenies(f))
														for f in families ]
	
	# 両親が補完されていない家系
	def extract_both_known_families(self) -> list[Family]:
		families = [ family for family in self.small_families
						if not family.is_self() and
							not self.is_imputed(family.mat) and
							not self.is_imputed(family.pat) and
							family.mat_known and family.pat_known and
							any(not self.is_imputed(prog)
											for prog in family.progenies) ]
		
		return [ Family(f.mat, f.pat, self.extract_unimputed_progenies(f))
														for f in families ]
	
	# 片親が補完されていて片親がunknownな家系
	def extract_one_imputed_families(self) -> list[KnownFamily]:
		families = [ family for family in self.small_families
							if not family.is_self() and
								((self.is_imputed(family.mat) and
													not family.pat_known) or
								 (self.is_imputed(family.pat) and
								  					not family.mat_known)) ]
		
		return [ KnownFamily(f.mat, f.pat, f.mat_known, f.pat_known, [prog])
							for f in families
							for prog in self.extract_unimputed_progenies(f) ]
	
	# 片親が補完されていなくて片親がunknownな家系
	def extract_one_known_parent_families(self) -> list[KnownFamily]:
		families = [ family for family in self.small_families
							if not family.is_self() and
								((family.mat_known and
					 						not self.is_imputed(family.mat) and
					 						not family.pat_known) or
								 (family.pat_known and
								 			not self.is_imputed(family.pat) and
								  			not family.mat_known)) and
						 		any(not self.is_imputed(prog)
											for prog in family.progenies) ]
		
		return [ KnownFamily(f.mat, f.pat, f.mat_known, f.mat_known,
									self.extract_unimputed_progenies(f))
														for f in families ]
	
	def extract_self_families(self) -> list[KnownFamily]:
		families = [ family for family in self.large_self_families
						if family.is_self() and
								any(not self.is_imputed(s)
									for s in family.samples()) ]
		
		return [ KnownFamily(f.mat, f.pat, f.mat_known, f.mat_known,
									f.progenies) for f in families ]
	
	def extract_small_self_families(self) -> list[KnownFamily]:
		families = [ family for family in self.small_families
						if family.is_self() and
								any(self.is_imputed(s)
									for s in family.samples()) and
								any(not self.is_imputed(s)
									for s in family.samples()) ]
		
		return [ KnownFamily(f.mat, f.pat, f.mat_known, f.mat_known,
									f.progenies) for f in families ]
	
	def extract_self_non_imputed_families(self) -> list[KnownFamily]:
		families = [ family for family in self.small_families
						if family.is_self() and
								all(not self.is_imputed(s)
									for s in family.samples()) ]
		
		return [ KnownFamily(f.mat, f.pat, f.mat_known, f.mat_known,
									f.progenies) for f in families ]
	
	def extract_self_parent_imputed_families(self) -> list[KnownFamily]:
		families = [ family for family in self.large_self_families
						if family.is_self() and self.is_imputed(family.mat) ]
		
		return [ KnownFamily(f.mat, f.pat, f.mat_known, f.mat_known,
									f.progenies) for f in families ]
	
	def extract_self_parent_non_imputed_families(self) -> list[KnownFamily]:
		families = [ family for family in self.large_self_families
						if family.is_self() and
							not self.is_imputed(family.mat) and
							any(not self.is_imputed(s)
									for s in family.samples()) ]
		
		return [ KnownFamily(f.mat, f.pat, f.mat_known, f.mat_known,
									f.progenies) for f in families ]
	
	# 両親は補完されていないが子どもの一部が補完されている家系
	def extract_progenies_imputed_families(self) -> list[KnownFamily]:
		families = [ family for family in self.small_families
								if not family.is_self() and
									(family.mat_known or family.pat_known) and
									not self.is_imputed(family.mat) and
									not self.is_imputed(family.pat) and
									any(self.is_imputed(prog)
											for prog in family.progenies) ]
		return [ KnownFamily(f.mat, f.pat, f.mat_known,
									f.mat_known, f.progenies)
														for f in families ]
	
	def extract_isolated_samples(self) -> list[str]:
		# 繋がっているサンプルがあっても、
		# 家系の全サンプルがphasingされていないなら孤立とみなす
		samples: set[str] = set()
		for family in self.small_families:
			for s in family.samples():
				if s != '0' and s not in self.imputed_samples:
					samples.add(s)
		return list(samples)
	
	def extract_non_imputed_samples(self) -> list[str]:
		samples = set(s for family in self.small_families
						for s in family.samples()
						if s != '0' and s not in self.imputed_samples)
		return list(samples)
	
	def display_info(self, out: TextIO) -> None:
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
	def make_families(ped: PedigreeTable, samples: list[str],
										lower_progs: int) -> list[KnownFamily]:
		families = ped.make_families(samples)
		new_families: list[KnownFamily] = []
		set_samples = set(samples)
		for family in families:
			# Families without any progeny in the VCF are deleted.
			filtered_progs: list[str] = []
			for prog in family.progenies:
				if prog in set_samples:
					filtered_progs.append(prog)
			if not filtered_progs:
				continue
			
			mat = family.mat
			pat = family.pat
			mat_known = mat in set_samples
			pat_known = pat in set_samples
			if len(filtered_progs) >= lower_progs and (mat_known or pat_known):
				new_families.append(KnownFamily(mat, pat, mat_known,
													pat_known, filtered_progs))
			else:
				# Treat parents who are not in the VCF as unknown
				mat_mod = mat if mat_known else '0'
				pat_mod = pat if pat_known else '0'
				new_families.append(KnownFamily(mat_mod, pat_mod, mat_known,
													pat_known, filtered_progs))
		
		return new_families
	
	@staticmethod
	def create(ped: PedigreeTable, samples: list[str], lower_progs: int,
								family_indices: list[int]) -> SampleManager:
		families = SampleManager.make_families(ped, samples, lower_progs)
		if family_indices:	# debug用にFamilyを絞って問題を小さくする
			families = [ families[i] for i in family_indices ]
		
		large_families = []
		large_self_families = []
		small_families = []
		for f in families:
			# 片親のGenotypeが無くても、後代の数が十分なら不明の親にしない
			if (not f.is_self() and f.num_progenies() >= lower_progs and
													not f.is_both_unknown()):
				large_families.append(f)
			elif (f.is_self() and f.num_progenies() >= lower_progs and
														f.is_both_known()):
				large_self_families.append(f)
			else:
				small_families.append(f)
		
		if not large_families and not large_self_families:
			raise SampleException()
		
		return SampleManager(ped, large_families, large_self_families,
													small_families, lower_progs)
