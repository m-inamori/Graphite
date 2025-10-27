
# coding: utf-8
# VCFFillableRecord.py

from __future__ import annotations
from functools import reduce
from itertools import *
from collections import Counter

import VCF
from VCFFamily import *
from GenoRecord import GenoRecord
from VCFHeteroHomo import VCFHeteroHomo
from VCFImpFamilyRecord import FillType, VCFImpFamilyRecord
from TypeDeterminer import ParentComb
from common import *
from Genotype import Genotype


#################### VCFFillableRecord ####################

class VCFFillableRecord(VCFFamilyRecord):
	def __init__(self, pos: int, geno: list[int], index: int,
									type: FillType, pair: ParentComb,
									probs: list[VCF.Probs]) -> None:
		super().__init__(pos, geno)
		self.index: int = index
		self.type: FillType = type
		self.pair: ParentComb = pair
		self.probs: list[VCF.Probs] = probs
	
	def gt_from_parent(self, mat_from: int, pat_from: int) -> int:
		if mat_from == 0 or pat_from == 0:
			return Genotype.NA
		
		gt_from_mat = self.get_mat_allele(mat_from-1)
		gt_from_pat = self.get_pat_allele(pat_from-1)
		return gt_from_mat | (gt_from_pat << 1) | 4
	
	# patからはどのアレルが来るのか分からない
	def gt_from_mat(self, mat_from: int, i: int) -> int:
		prog_gt = self.unphased(i)
		# matからはこのアレルが来ているはず
		mat_allele = (self.mat_gt() >> (mat_from-1)) & 1
		pat_gt = self.pat_gt()
		if prog_gt == Genotype.UN_00:
			if mat_allele == 0:
				return Genotype.PH_00
			elif pat_gt < Genotype.PH_11:
				return Genotype.PH_10
			else:
				return Genotype.PH_11
		elif prog_gt == Genotype.UN_11:
			if mat_allele == 1:
				return Genotype.PH_11
			elif pat_gt == Genotype.PH_00:
				return Genotype.PH_00
			else:
				return Genotype.PH_01
		elif prog_gt == Genotype.UN_01:
			if mat_allele == 0:
				return Genotype.PH_01
			else:
				return Genotype.PH_10
		else:
			return Genotype.NA
	
	def gt_from_pat(self, pat_from: int, i: int) -> int:
		prog_gt = self.unphased(i)
		# matからはこのアレルが来ているはず
		pat_allele = (self.pat_gt() >> (pat_from-1)) & 1
		mat_gt = self.mat_gt()
		if prog_gt == Genotype.UN_00:
			if pat_allele == 0:
				return Genotype.PH_00
			elif mat_gt == Genotype.PH_11:
				return Genotype.PH_11
			else:
				return Genotype.PH_01
		elif prog_gt == Genotype.UN_11:
			if pat_allele == 1:
				return Genotype.PH_11
			elif mat_gt == Genotype.PH_00:
				return Genotype.PH_10
			else:
				return Genotype.PH_00
		elif prog_gt == Genotype.UN_01:
			if pat_allele == 0:
				return Genotype.PH_10
			else:
				return Genotype.PH_01
		else:
			return Genotype.NA
	
	def mat_from(self, i: int) -> int:
		if self.is_NA(i):
			return 0
		elif not self.is_mat_hetero():
			return 0
		else:
			a = self.get_mat_allele(0)
			# 文字列にしたときに、1文字目がどうなるか
			gt = self.geno[i]
			if gt == Genotype.UN_11 or gt == Genotype.PH_10 or gt == Genotype.PH_11:
				b = 1
			else:
				b = 0
			
			if a == b:
				return 1
			else:
				return 2
	
	def pat_from(self, i: int) -> int:
		if self.is_NA(i):
			return 0
		elif not self.is_pat_hetero():
			return 0
		else:
			# ここおかしい
			a = self.get_pat_allele(0)
			# 文字列にしたときに、3文字目がどうなるか
			gt = self.geno[i]
			if gt == Genotype.UN_00 or gt == Genotype.PH_00 or gt == Genotype.PH_10:
				b = 0
			else:
				b = 1
			
			if a == b:
				return 1
			else:
				return 2
	
	def inverse_prog_gt(self, gt: int, inv_mat: bool, inv_pat: bool) -> int:
		def inverse_allele(a: int, inv: bool) -> int:
			return 1 - a if inv else a
		
		allele1 = Genotype.get_allele(gt, 0)
		allele2 = Genotype.get_allele(gt, 1)
		return (4 | inverse_allele(allele1, inv_mat) |
					(inverse_allele(allele2, inv_pat) << 1))
	
	def inverse_prog_gts(self, prog_gts: list[int],
							inv_mat: bool, inv_pat: bool) -> list[int]:
		inv_prog_gts = [ self.inverse_prog_gt(gt, inv_mat, inv_pat)
												for gt in prog_gts ]
		return inv_prog_gts
	
	def inverse_parents_gts(self, inv_mat: bool, inv_pat: bool) -> None:
		# both must be hetero
		if inv_mat:
			self.geno[0] = Genotype.inverse(self.geno[0])
		if inv_pat:
			self.geno[1] = Genotype.inverse(self.geno[1])
	
	def is_near_prog_gts(self, gts: list[int]) -> bool:
		num = 0
		dist = 0
		for i in range(len(gts)):
			if self.geno[i+2] != Genotype.UN_01:
				num += 1
			# self.genoはnon-phasedの前提だがそれでよいか？
			if Genotype.unphased(gts[i]) != self.geno[i+2]:
				dist += 1
		return dist < max(1, num // 2)
	
	def modify_gts(self, new_prog_gts: list[int]) -> None:
		for inv_mat, inv_pat in product((False, True), repeat=2):
			inv_prog_gts = self.inverse_prog_gts(new_prog_gts, inv_mat, inv_pat)
			if self.is_near_prog_gts(inv_prog_gts):
				self.inverse_parents_gts(inv_mat, inv_pat)
				for i in range(2, len(self.geno)):
					self.geno[i] = inv_prog_gts[i-2]
				return
		for i in range(2, len(self.geno)):
			self.geno[i] = new_prog_gts[i-2]
	
	def modify_parents_type(self) -> None:
		mat_gt = self.mat_gt()
		pat_gt = self.pat_gt()
		if (self.pair != ParentComb.P00x11 and
				((mat_gt == Genotype.PH_00 and pat_gt == Genotype.PH_11) or
				 (mat_gt == Genotype.PH_11 and pat_gt == Genotype.PH_00))):
			self.pair = ParentComb.P00x11
	
	def from_which_chrom(self, i: int, mat: bool) -> int:
		j = 0 if mat else 1
		parent_allele = self.get_allele(j, 0)
		return 1 if parent_allele == self.get_allele(i, j) else 2
	
	@staticmethod
	def convert(records: list[VCFImpFamilyRecord],
								vcf: VCFHeteroHomo) -> list[VCFFillableRecord]:
		ref_vcf = vcf.vcf
		cols = ref_vcf.extract_columns(vcf.samples)
		new_records: list[VCFFillableRecord] = []
		for i, record in enumerate(records):
			type = record.get_fill_type()
			probs = ref_vcf.records[i].parse_PL(record.geno, cols)
			r = VCFFillableRecord(record.pos, record.geno,
										record.index, type, record.pair, probs)
			new_records.append(r)
		return new_records
	
	@staticmethod
	def merge(records: list[VCFFillableRecord],
								samples: list[str]) -> GenoRecord:
		pos = records[0].pos
		geno = records[0].geno[:]
		for record in records[1:]:
			geno.extend(record.geno)
		return GenoRecord(pos, geno)
	
	def decide_by_majority(self, GTs: list[int]) -> int:
		# 多数決できるなら多数決
		c = Counter(GTs)
		dic_num = Counter(c.values())
		max_num = max(dic_num.keys())
		if dic_num[max_num] == 1:	# 最多のGenotypeは一つ
			for GT, num in c.items():
				if num == max_num:
					return GT
			else:
				# ここには来ないはず
				return GT
		else:						# 複数ある場合は乱数的なものを使う
			candidate_GTs = []
			for GT, num in c.items():
				if num == max_num:
					candidate_GTs.append(GT)
			candidate_GTs.sort()
			# posをhash化してGTを選ぶ
			d = len(candidate_GTs)
			s = 0
			pos = self.pos
			while pos > 0:
				s += pos % 10
				pos //= 10
			i = s % d
			return candidate_GTs[i]
	
	# 親が00x11のときGTと違うなら0|0と1|1を入れ替える
	def swap_parents(self, pos: int, GT: int) -> None:
		if GT == self.geno[pos]:
			return
		elif GT in (Genotype.PH_00, Genotype.PH_11):
			is_mat_00 = (pos == 0) ^ (GT == Genotype.PH_11)
			self.geno[0] = Genotype.PH_00 if is_mat_00 else Genotype.PH_11
			self.geno[1] = Genotype.PH_11 if is_mat_00 else Genotype.PH_00
			# 子どもも入れ換えなければならない
			prog_GT = Genotype.PH_01 if is_mat_00 else Genotype.PH_10
			for i in range(2, len(self.geno)):
				self.geno[i] = prog_GT
	
	@staticmethod
	def decide_duplicated_Genotype(records: list[VCFFillableRecord],
											ps: list[tuple[int, int]]) -> int:
		# 全部./.ならそのまま
		# 子どもがいたらそれ
		# ./.を除いて全部同じだったらそれ
		# ./.と0/0 x 1/1で全部なら0/0 x 1/1で多数決
		# ./.と0/0 x 1/1を除いて全部同じだったらそれ
		# ./.と0/0 x 1/1を除いて複数種あれば多数決
		GTs = [ records[i].geno[j] for i, j in ps ]
		if all(GT == Genotype.NA for GT in GTs):
			return Genotype.NA
		
		# 子どもが優先
		for (i, j), GT in zip(ps, GTs):
			if j >= 2 and GT != Genotype.NA:
				return GT
		
		GTs2 = [ GT for GT in GTs if GT != Genotype.NA ]
		if is_all_same(GTs2):
			return GTs2[0]
		
		GTs3 = [ GT for (i, j), GT in zip(ps, GTs)
						if GT != Genotype.NA and
						   records[i].pair != ParentComb.P00x11 ]
		if not GTs3:	# 0/0 x 1/1しかない
			return records[0].decide_by_majority(GTs2)
		elif is_all_same(GTs3):
			return GTs3[0]
		else:
			return records[0].decide_by_majority(GTs3)
	
	@staticmethod
	def integrate_each_sample(records: list[VCFFillableRecord],
										ps: list[tuple[int, int]]) -> None:
		GT = VCFFillableRecord.decide_duplicated_Genotype(records, ps)
		
		for i, j in ps:
			# 0/0 x 1/1の親なら相手と入れ替えなければならない
			if records[i].pair == ParentComb.P00x11 and j <= 1:
				records[i].swap_parents(j, GT)
	
	@staticmethod
	def integrate(records: list[VCFFillableRecord], samples: list[str],
						pos_samples: list[list[tuple[int, int]]]) -> GenoRecord:
		pos = records[0].pos
		# 0|0 x 1|1ならひっくり返せる
		for ps in pos_samples:
			GTs = [ records[i].geno[j] for i, j in ps ]
			if not is_all_same(GTs):
				VCFFillableRecord.integrate_each_sample(records, ps)
		
		# ひっくり返した後にGTをまとめる
		geno: list[int] = []
		for ps in pos_samples:
			i, j = ps[0]
			geno.append(records[i].geno[j])
		return GenoRecord(pos, geno)
