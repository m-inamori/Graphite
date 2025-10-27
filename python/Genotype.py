from __future__ import annotations
from itertools import product

# coding: utf-8
# Genotype.py


#################### Genotype ####################

class Genotype:
	UN_00 = 0	# 0/0
	UN_01 = 1	# 0/1
	UN_11 = 2	# 1/1
	NA = 3
	PH_00 = 4	# 0|0
	PH_10 = 5	# 1|0
	PH_01 = 6	# 0|1
	PH_11 = 7	# 1|1
	
	def __init__(self, s: str):
		self.phasing: bool	= s[1] == '|'
		self.gt1: str		= s[0]
		self.gt2: str		= s[2]
	
	def gts(self) -> list[str]:
		return [self.gt1, self.gt2]
	
	def __str__(self) -> str:
		return self.gt1 + ('|' if self.phasing else '/') + self.gt2
	
	@staticmethod
	def is_ref_homo(gt: int) -> bool:
		return gt == Genotype.UN_00 or gt == Genotype.PH_00
	
	@staticmethod
	def is_alt_homo(gt: int) -> bool:
		return gt == Genotype.UN_11 or gt == Genotype.PH_11
	
	@staticmethod
	def is_hetero(gt: int) -> bool:
		return (gt == Genotype.UN_01 or gt == Genotype.PH_10 or
										gt == Genotype.PH_01)
	
	@staticmethod
	def unphased(gt: int) -> int:
		if gt < 4:
			return gt
		elif gt == 4:
			return 0
		elif gt < 7:
			return 1
		else:
			return 2
	
	@staticmethod
	def int_to_gt(n: int) -> str:
		if n == 0:
			return '0/0'
		elif n == 1:
			return '0/1'
		elif n == 2:
			return '1/1'
		else:
			return './.'
	
	@staticmethod
	def conflicts(mat: Genotype, pat: Genotype, prog: Genotype,
									considers_phasing: bool = True) -> bool:
		if considers_phasing and prog.phasing:
			return not (prog.gt1 in mat.gts() and prog.gt2 in pat.gts())
		else:
			if prog.gt1 == prog.gt2:
				return not (prog.gt1 in mat.gts() and prog.gt2 in pat.gts())
			else:
				return not ((prog.gt1 in mat.gts() and prog.gt2 in pat.gts()) or
							(prog.gt1 in pat.gts() and prog.gt2 in mat.gts()))
	
	@staticmethod
	def is_valid(gt: int, mat_gt: int, pat_gt: int) -> bool:
		mat_gts = Genotype.possible_gts(mat_gt)
		pat_gts = Genotype.possible_gts(pat_gt)
		for mat_gt, pat_gt in product(mat_gts, pat_gts):
			if mat_gt + pat_gt == 0 and Genotype.is_00(gt):
				return True
			elif mat_gt + pat_gt == 1 and Genotype.is_01(gt):
				return True
			elif mat_gt + pat_gt == 2 and Genotype.is_11(gt):
				return True
		return False
	
	@staticmethod
	def possible_gts(gt: int) -> list[int]:
		if Genotype.is_00(gt):
			return [0]
		elif Genotype.is_11(gt):
			return [1]
		else:
			return [0, 1]
	
	@staticmethod
	def gt_to_int(gt: str) -> int:
		if '.' in gt[:3]:
			return 3
		gt1 = 0 if gt[0] == '0' else 1
		gt2 = 0 if gt[2] == '0' else 1
		return gt1 + gt2
	
	@staticmethod
	def phased_gt_to_int(gt: str) -> int:
		gt1 = 0 if gt[0] == '0' else 1
		gt2 = 0 if gt[2] == '0' else 1
		return gt1 | (gt2 << 1)
	
	# N/A: 3
	# non-phased: 0-2
	# phased: 4-7
	@staticmethod
	def all_gt_to_int(gt: str) -> int:
		if '.' in gt[:3]:
			return 3
		elif gt[1] == '/':
			return Genotype.gt_to_int(gt)
		else:
			return Genotype.phased_gt_to_int(gt) | 4
	
	# 0-7 -> 0-3
	@staticmethod
	def all_int_gt_to_int_gt(gt_int: int) -> int:
		if gt_int < 5:
			return gt_int & 3
		else:
			# 6(0|1) -> 1
			# 7(1|1) -> 2
			return gt_int - 5
	
	# a1, a2 -> phased genotype
	@staticmethod
	def from_alleles(a1: int, a2: int) -> int:
		return a1 | (a2 << 1) | 4
	
	@staticmethod
	def is_NA(gt_int: int) -> bool:
		return gt_int == Genotype.NA
	
	@staticmethod
	def is_00(gt_int: int) -> bool:
		return (gt_int & 3) == 0
	
	@staticmethod
	def is_01(gt_int: int) -> bool:
		return (gt_int == Genotype.UN_01 or
				gt_int == Genotype.PH_01 or gt_int == Genotype.PH_10)
	
	@staticmethod
	def is_11(gt_int: int) -> bool:
		return gt_int == Genotype.UN_11 or gt_int == Genotype.PH_11
	
	@staticmethod
	def is_phased(gt_int: int) -> bool:
		return gt_int >= 4
	
	# When de-phased, can they be considered to be the same genotype?
	@staticmethod
	def is_same_non_phased_genotype(gt1: int, gt2: int) -> bool:
		return (Genotype.all_int_gt_to_int_gt(gt1) == 
				Genotype.all_int_gt_to_int_gt(gt2))
	
	@staticmethod
	def int_to_phased_gt(gt_int: int) -> str:
		if gt_int == 0:
			return '0|0'
		elif gt_int == 1:
			return '1|0'
		elif gt_int == 2:
			return '0|1'
		else:
			return '1|1'
	
	@staticmethod
	def int_to_all_gt(gt_int: int) -> str:
		if gt_int < 4:
			return Genotype.int_to_gt(gt_int)
		else:
			return Genotype.int_to_phased_gt(gt_int & 3)
	
	@staticmethod
	def is_homo(all_gt_int: int) -> bool:
		if all_gt_int < 6:
			return (all_gt_int & 1) == 0
		else:
			return all_gt_int == 7
	
	# 0|1 <-> 1|0 phased前提
	@staticmethod
	def inverse(gt: int) -> int:
		if gt == Genotype.PH_01:
			return Genotype.PH_10
		elif gt == Genotype.PH_10:
			return Genotype.PH_01
		else:
			return gt
	
	# phasedのGenotypeのj(0 or 1)側のallele
	@staticmethod
	def get_allele(gt: int, j: int) -> int:
		return (gt >> j) & 1
	
	# ハプロタイプを決めたときのGenotype
	@staticmethod
	def gt_by_haplotypes(hc1: int, hc2: int,
							mat_gt: int, pat_gt: int) -> int:
		return ((mat_gt >> hc1) & 1) | (((pat_gt >> hc2) & 1) << 1)
	
	@staticmethod
	def find_key_position(info: str, key: str) -> int:
		for i, t in enumerate(info.split(':')):
			if t == key:
				return i
		else:
			return -1
