from __future__ import annotations

# coding: utf-8
# VCFFillable.py

from itertools import *
from collections import defaultdict

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from VCFFamily import *
from VCFImpFamilyRecord import FillType, VCFImpFamilyRecord
from VCFFillableRecord import VCFFillableRecord
from VCFHeteroHomo import *
from group import Groups
from RecordSet import RecordSet
from Genotype import Genotype
from option import *
from common import *


#################### VCFFillable ####################

class VCFFillable(VCFFamilyBase):
	def __init__(self, samples: list[str], records: list[VCFFillableRecord],
														vcf: VCFSmall) -> None:
		VCFFamilyBase.__init__(self, samples, vcf)
		self.records: list[VCFFillableRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def modify(self, is_phased_changable: bool) -> None:
		# FillTypeでrecordを分ける
		groups = Groups.create(self.records)
		for record_set in groups.generate_record_sets():
			record_set.determine_parents_phasing()
			self.impute_core(record_set)
		
		for i, record in enumerate(self.records):
			# 家系ごとで./.にしたGenotypeを補完
			if record.type == FillType.MAT:
				self.impute_NA_mat(i)
			elif record.type == FillType.PAT:
				self.impute_NA_pat(i)
			elif record.type in (FillType.IMPUTABLE, FillType.UNABLE):
				self.impute_others(i)
	
	def phase_hetero_hetero(self) -> None:
		# FillTypeでrecordを分ける
		groups = Groups.create(self.records)
		for record_set in groups.generate_record_sets():
			record_set.impute(False)
		
		# この部分はフルで必要？
		for i, record in enumerate(self.records):
			# 家系ごとで./.にしたGenotypeを補完
			if record.type == FillType.MAT:
				self.impute_NA_mat(i)
			elif record.type == FillType.PAT:
				self.impute_NA_pat(i)
			elif record.type in (FillType.IMPUTABLE, FillType.UNABLE):
				self.impute_others(i)
	
	def impute_core(self, record_set: RecordSet) -> None:
		def gts(record: Optional[VCFFillableRecord]) -> list[int]:
			N = len(self.samples) - 2
			return record.geno[2:] if record else [Genotype.NA] * N
		
		def both_nearests() -> tuple[int, int]:
			# ここに来るときは、全てNoneではない
			if (record is None or mat_record1 is None or mat_record2 is None
							   or pat_record1 is None or pat_record2 is None):
				return (0, 0)
			if record_set.is_prev_nearer(True):
				mat_from = prev_mat_from
			else:
				mat_from = next_mat_from
			if record.pos * 2 < pat_record1.pos + pat_record2.pos:
				pat_from = prev_pat_from
			else:
				pat_from = next_pat_from
			return (mat_from, pat_from)
		
		# [(mat_from, pat_from)] -> (mat_from, pat_from)
		def nearest_froms(pairs: list[tuple[int, int]]) -> tuple[int, int]:
			if len(pairs) == 4:		# N/Aのときにまれにあり得る
				return both_nearests()
			if pairs[0][0] == pairs[1][0]:		# matが同じ
				if record is None or pat_record1 is None or pat_record2 is None:
					return (0, 0)
				if record.pos * 2 < pat_record1.pos + pat_record2.pos:
					return (pairs[0][0], prev_pat_from)
				else:
					return (pairs[0][0], next_pat_from)
			elif pairs[0][1] == pairs[1][1]:	# patが同じ
				if record is None or mat_record1 is None or mat_record2 is None:
					return (0, 0)
				if record.pos * 2 < mat_record1.pos + mat_record2.pos:
					return (prev_mat_from, pairs[0][1])
				else:
					return (next_mat_from, pairs[0][1])
			else:	# 両親とも乗り換えている（滅多にない）
				return both_nearests()
		
		def select_pair(pairs: list[tuple[int, int]], gt: int,
									selected: bool = False) -> tuple[int, int]:
			if record is None:
				return (0, 0)
			elif not pairs:
				return (0, 0)
			elif len(pairs) == 1:
				return pairs[0]
			elif gt == Genotype.NA:
				return nearest_froms(pairs)
			elif selected:
				return nearest_froms(pairs)
			else:
				new_pairs = [ v for v in pairs
							if Genotype.unphased(gt) ==
								Genotype.unphased(record.gt_from_parent(*v)) ]
				pair = select_pair(new_pairs, gt, True)
				if pair != (0, 0):
					return pair
				else:
					return select_pair(pairs, gt, True)
		
		def is_same_gts(gt1: int, gt2: int) -> bool:
			if gt2 == Genotype.UN_01:
				return gt1 == Genotype.PH_01 or gt1 == Genotype.PH_10
			elif gt2 == Genotype.UN_00:
				return gt1 == Genotype.PH_00
			elif gt2 == Genotype.UN_11:
				return gt1 == Genotype.PH_11
			else:
				return False
		
		def is_near(gts1: list[int], gts2: list[int]) -> bool:
			num = sum(1 for gt2 in gts2 if gt2 != Genotype.UN_01)
			dist = sum(1 for gt1, gt2 in zip(gts1, gts2)
											if not is_same_gts(gt1, gt2))
			return dist < max(1, num // 2)
		
		def distance(gts1: list[int], gts2: list[int]) -> int:
			return sum(1 for gt1, gt2 in zip(gts1, gts2)
											if not is_same_gts(gt1, gt2))
		
		def inverse_gt(gt: int, inv_mat: bool, inv_pat: bool) -> int:
			def inverse(a: int) -> int:
				return 1 if a == 0 else 0
			
			# phasingされているはずなので、ここは適当
			if Genotype.is_NA(gt):
				return Genotype.NA
			
			if gt == 0:
				gt = 4
			elif gt == 1:
				gt = 6
			elif gt == 2:
				gt = 7
			
			a1 = Genotype.get_allele(gt, 0)
			a2 = Genotype.get_allele(gt, 1)
			if inv_mat:
				a1 = inverse(a1)
			if inv_pat:
				a2 = inverse(a2)
			return Genotype.from_alleles(a1, a2)
		
		def inverse_gts(gts: list[int], inv_mat: bool,
										inv_pat: bool) -> list[int]:
			return [ inverse_gt(gt, inv_mat, inv_pat) for gt in gts ]
		
		def inverse_parents_gts(record: VCFFillableRecord,
								inv_mat: bool, inv_pat: bool) -> None:
			# both must be hetero
			if inv_mat:
				record.geno[0] = Genotype.inverse(record.geno[0])
			if inv_pat:
				record.geno[1] = Genotype.inverse(record.geno[1])
		
		def modify_gts(gts: list[int],
						record: Optional[VCFFillableRecord]) -> None:
			if record is None:
				return
			
			orig_gts = record.geno[2:]
			# 元のGenotypeに最も近いinverseの組合せを使う
			min_dist = distance(gts, orig_gts)
			min_gts = gts
			min_i = 0
			# gtsがnon-phasedのときがあるから、適当にphasingしてinverseする
			for i in range(4):
				inv_mat = (i & 2) == 2
				inv_pat = (i & 1) == 1
				new_gts = inverse_gts(gts, inv_mat, inv_pat)
				dist = distance(new_gts, orig_gts)
				if dist < min_dist:
					min_dist = dist
					min_gts = new_gts
					min_i = i
			
			inv_mat = (min_i & 2) == 2
			inv_pat = (min_i & 1) == 1
			inverse_parents_gts(record, inv_mat, inv_pat)
			for i in range(2, len(record.geno)):
				record.geno[i] = min_gts[i-2]
			
			"""
			for inv_mat, inv_pat in product((False, True), repeat=2):
				new_gts = inverse_gts(gts, inv_mat, inv_pat)
				if is_near(new_gts, orig_gts):
					inverse_parents_gts(record, inv_mat, inv_pat)
					for i in range(2, len(record.geno)):
						record.geno[i] = new_gts[i-2]
					break
			else:
				new_gts = gts
			for i in range(2, len(record.geno)):
				record.geno[i] = new_gts[i-2]
			"""
		
		# どちらから来たか決める
		def select_from(froms: list[int], record1: Optional[GenoRecord],
										  record2: Optional[GenoRecord]) -> int:
			if len(froms) == 1:
				return froms[0]
			
			# どちらかは0でない前提
			if froms[0] == 0:
				return froms[1]
			elif froms[1] == 0:
				return froms[0]
			else:	# 両側にRecordがある
				# 近い方を選ぶ
				if record is None or record1 is None or record2 is None:
					# ここには来ないはず
					return 0
				elif record.pos * 2 < record1.pos + record2.pos:
					return froms[0]
				else:
					return froms[1]
		
		def modify_parents_type(record: VCFFillableRecord) -> None:
			if record is not None:
				record.modify_parents_type()
		
		record, mat_record1, mat_record2, pat_record1, pat_record2 = \
															record_set.records()
		if record is None:
			return
		new_gts = record.geno[2:]
		for i, gt, mat_gt1, mat_gt2, pat_gt1, pat_gt2 in zip(count(2),
											*map(gts, record_set.records())):
			prev_mat_from = record_set.from_which_chrom_prev_mat(mat_gt1)
			next_mat_from = record_set.from_which_chrom_next_mat(mat_gt2)
			prev_pat_from = record_set.from_which_chrom_prev_pat(pat_gt1)
			next_pat_from = record_set.from_which_chrom_next_pat(pat_gt2)
			if ((prev_mat_from == 0 and next_mat_from == 0) or
					(prev_pat_from == 0 and next_pat_from == 0)):
				new_gts[i-2] = record.geno[i]
				continue
			
			mat_froms = unique_list(prev_mat_from, next_mat_from)
			pat_froms = unique_list(prev_pat_from, next_pat_from)
			pairs = [ x for x in product(mat_froms, pat_froms)
										if x[0] != 0 and x[1] != 0 ]
			mat_from, pat_from = select_pair(pairs, gt)
			if (mat_from, pat_from) == (0, 0):
				if any(mat_from != 0 for mat_from in mat_froms):
					mat_from_nz = select_from(mat_froms,
												mat_record1, mat_record2)
					new_gts[i-2] = record.gt_from_mat(mat_from_nz, i)
				elif any(pat_from != 0 for pat_from in pat_froms):
					pat_from_nz = select_from(pat_froms,
												pat_record1, pat_record2)
					new_gts[i-2] = record.gt_from_pat(pat_from_nz, i)
				else:
					continue	# phasingしない
			else:
				new_gts[i-2] = record.gt_from_parent(mat_from, pat_from)
		
		modify_gts(new_gts, record)		# 両親が間違っていた時の処理
		modify_parents_type(record)
	
	def find_prev_same_type_record(self, i: int, k: int
											) -> Optional[VCFFillableRecord]:
		type = self.records[i].type
		for j in range(i - 1, -1, -1):
			record = self.records[j]
			if record.type == type and not record.is_NA(k):
				return record
		else:
			return None
	
	def find_next_same_type_record(self, i: int, k: int
											) -> Optional[VCFFillableRecord]:
		type = self.records[i].type
		for j in range(i + 1, len(self.records)):
			record = self.records[j]
			if record.type == type and not record.is_NA(k):
				return record
		else:
			return None
	
	# 家系ごとで./.にしたGenotypeを補完
	def impute_NA_mat(self, i: int) -> None:
		record = self.records[i]
		for k in range(2, len(record.geno)):
			if record.is_NA(k):
				self.impute_NA_mat_each(i, k)
	
	def create_recordset_mat(self, record: Optional[VCFFillableRecord],
						prev_record: Optional[VCFFillableRecord],
						next_record: Optional[VCFFillableRecord]) -> RecordSet:
		return RecordSet(record, prev_record, next_record, None, None)
	
	def create_recordset_pat(self, record: Optional[VCFFillableRecord],
						prev_record: Optional[VCFFillableRecord],
						next_record: Optional[VCFFillableRecord]) -> RecordSet:
		return RecordSet(record, None, None, prev_record, next_record)
	
	def impute_NA_mat_each(self, i: int, k: int) -> None:
		def select_mat(pairs: list[tuple[int, int]]) -> int:
			if len(pairs) == 1:
				return pairs[0][0]
			elif prev_record is None and next_record is None:
				return 0
			elif next_record is None:
				return prev_mat_from
			elif prev_record is None:
				return next_mat_from
			elif record.pos * 2 < prev_record.pos + next_record.pos:
				return prev_mat_from
			else:
				return next_mat_from
		
		record = self.records[i]
		prev_record = self.find_prev_same_type_record(i, k)
		next_record = self.find_next_same_type_record(i, k)
		recordset = self.create_recordset_mat(record, prev_record, next_record)
		prev_gt = (prev_record.geno[k] if prev_record is not None
										else Genotype.NA)
		next_gt = (next_record.geno[k] if next_record is not None
										else Genotype.NA)
		prev_mat_from = recordset.from_which_chrom_prev_mat(prev_gt)
		next_mat_from = recordset.from_which_chrom_next_mat(next_gt)
		mat_froms = unique_list(prev_mat_from, next_mat_from)
		pairs = [ (x, 1) for x in mat_froms if x != 0 ]
		if not pairs:	# 両側ともにmatが無い（まず無い）
			return
		
		mat_from = select_mat(pairs)
		if mat_from == 0:
			return
		
		record.geno[k] = record.gt_from_parent(mat_from, 1)
	
	# 家系ごとで./.にしたGenotypeを補完
	def impute_NA_pat(self, i: int) -> None:
		record = self.records[i]
		for k in range(2, len(record.geno)):
			if record.is_NA(k):
				self.impute_NA_pat_each(i, k)
	
	def impute_NA_pat_each(self, i: int, k: int) -> None:
		def select_pat(pairs: list[tuple[int, int]]) -> int:
			if len(pairs) == 1:
				return pairs[0][1]
			elif prev_record is None and next_record is None:
				return 0
			elif next_record is None:
				return prev_pat_from
			elif prev_record is None:
				return next_pat_from
			elif record.pos * 2 < prev_record.pos + next_record.pos:
				return prev_pat_from
			else:
				return next_pat_from
		
		record = self.records[i]
		prev_record = self.find_prev_same_type_record(i, k)
		next_record = self.find_next_same_type_record(i, k)
		recordset = self.create_recordset_pat(record, prev_record, next_record)
		prev_gt = (prev_record.geno[k]
					if prev_record is not None else Genotype.NA)
		next_gt = (next_record.geno[k]
					if next_record is not None else Genotype.NA)
		prev_pat_from = recordset.from_which_chrom_prev_pat(prev_gt)
		next_pat_from = recordset.from_which_chrom_next_pat(next_gt)
		pat_froms = unique_list(prev_pat_from, next_pat_from)
		pairs = [ (1, x) for x in pat_froms if x != 0 ]
		if not pairs:	# 両側ともにpatが無い（まず無い）
			return
		
		pat_from = select_pat(pairs)
		record.geno[k] = record.gt_from_parent(1, pat_from)
	
	def impute_others(self, i: int) -> None:
		# GTはphased
		def correct(GT: int, i: int, record: VCFFillableRecord) -> int:
			if record.is_NA(i):
				return (record.pos % 2) | ((record.pos >> 1) % 2) | 4
			elif record.is_00(i):
				if (GT & 1) == 1:		# 母側が1
					return Genotype.PH_10
				elif (GT & 2) == 2:		# 父側が1
					return Genotype.PH_01
				else:
					return Genotype.PH_00
			elif record.is_01(i):
				if (GT & 1) == 1:		# 母側が1
					return Genotype.PH_10
				elif (GT & 2) == 2:		# 父側が1
					return Genotype.PH_01
				else:
					# どちらにすべきか判別がつかないので適当に選ぶ
					return (Genotype.PH_01 if record.pos % 2 == 0
												else Genotype.PH_10)
			else:
				if (GT & 1) == 0:		# 母側が0
					return Genotype.PH_01
				elif (GT & 2) == 0:		# 父側が0
					return Genotype.PH_10
				else:
					return Genotype.PH_11
		
		record = self.records[i]
		ks = [ k for k in range(2, len(record.geno))
								if not record.is_phased(k) ]
		if not ks:
			return
		
		mat_homo = record.is_homo(0)
		pat_homo = record.is_homo(1)
		for k in ks:
			mat_from = 1 if mat_homo else self.__find_mat_from(i, k)
			pat_from = 1 if pat_homo else self.__find_pat_from(i, k)
			GT = record.gt_from_parent(mat_from, pat_from)
			if Genotype.is_NA(GT):
				GT = correct(GT, k, record)
			record.geno[k] = GT
	
	def __find_mat_from(self, i: int, k: int) -> int:
		return self.__select_from(self.__find_prev_mat_from(i, k),
									self.__find_next_mat_from(i, k), i)
	
	def __find_pat_from(self, i: int, k: int) -> int:
		return self.__select_from(self.__find_prev_pat_from(i, k),
									self.__find_next_pat_from(i, k), i)
	
	def __select_from(self, f1: tuple[int, int],
							f2: tuple[int, int], i: int) -> int:
		i1, from1 = f1
		i2, from2 = f2
		if from1 == 0 and from2 == 0:
			# 前後がないとき乱数的に決める
			r0 = self.records[i]
			return r0.pos % 2 + 1
		
		if from1 == from2:
			return from1
		elif from2 == 0:
			return from1
		elif from1 == 0:
			return from2
		else:
			# 最後は物理距離で決める
			r0 = self.records[i]
			r1 = self.records[i1]
			r2 = self.records[i2]
			if r0.pos * 2 <= r2.pos + r1.pos:
				return from1
			else:
				return from2
	
	def __find_prev_mat_from(self, i: int, k: int) -> tuple[int, int]:
		for j in range(i - 1, -1, -1):
			from1 = self.records[j].mat_from(k)
			if from1 != 0:
				return (j, from1)
		return (-1, 0)
	
	def __find_next_mat_from(self, i: int, k: int) -> tuple[int, int]:
		for j in range(i + 1, len(self.records)):
			from1 = self.records[j].mat_from(k)
			if from1 != 0:
				return (j, from1)
		return (-1, 0)
	
	def __find_prev_pat_from(self, i: int, k: int) -> tuple[int, int]:
		for j in range(i - 1, -1, -1):
			from1 = self.records[j].pat_from(k)
			if from1 != 0:
				return (j, from1)
		return (-1, 0)
	
	def __find_next_pat_from(self, i: int, k: int) -> tuple[int, int]:
		for j in range(i + 1, len(self.records)):
			from1 = self.records[j].pat_from(k)
			if from1 != 0:
				return (j, from1)
		return (-1, 0)
	
	def remove_parents(self) -> VCFGeno:
		samples: list[str] = self.samples[2:]
		records: list[GenoRecord] = []
		for record in self.records:
			geno = record.geno[2:]
			r = GenoRecord(record.pos, geno)
			records.append(r)
		return VCFGeno(samples, records, self.vcf)
	
	@staticmethod
	def fill(vcfs: list[VCFHeteroHomo],
					records: list[VCFImpFamilyRecord]) -> VCFFillable:
		merged_records = VCFFillable.merge_records(vcfs, records)
		vcf: VCFFillable = VCFFillable(vcfs[0].samples, merged_records,
															vcfs[0].vcf)
		vcf.modify(True)
		return vcf
	
	@staticmethod
	def merge_records(vcfs: list[VCFHeteroHomo],
						records: list[VCFImpFamilyRecord],
						) -> list[VCFFillableRecord]:
		all_records: list[VCFImpFamilyRecord] = []
		for vcf in vcfs:
			for r in vcf.records:
				all_records.append(r)
		all_records.extend(records)
		all_records.sort(key=lambda record: record.index)
		return VCFFillableRecord.convert(all_records, vcfs[0])
	
	@staticmethod
	def collect_records(vcfs: list[VCFFillable]
									) -> list[list[VCFFillableRecord]]:
		rss: list[list[VCFFillableRecord]] = []
		if not vcfs:
			return rss
		
		for i in range(len(vcfs[0])):
			rss.append([ vcf.records[i] for vcf in vcfs ])
		return rss
	
	@staticmethod
	def integrate_samples(sss: list[list[str]], orig_samples: list[str]
							) -> tuple[list[str], list[list[tuple[int, int]]]]:
		dic = defaultdict(list)
		for i, samples in enumerate(sss):
			for j, sample in enumerate(samples):
				dic[sample].append((i, j))
		samples = [ sample for sample in orig_samples if sample in dic ]
		pos_samples = [ dic[sample] for sample in samples ]
		return (samples, pos_samples)
	
	# 重複したサンプルが一つになるようにVCFを統合する
	@staticmethod
	def integrate(ref_vcf: VCFSmall, rss: list[list[VCFFillableRecord]],
					sss: list[list[str]], orig_samples: list[str]) -> VCFGeno:
		samples, pos_samples = VCFFillable.integrate_samples(sss, orig_samples)
		records: list[GenoRecord] = []
		for rs in rss:
			records.append(VCFFillableRecord.integrate(rs, samples, pos_samples))
		return VCFGeno(samples, records, ref_vcf)
	
	@staticmethod
	def merge(vcfs: list[VCFFillable], orig_samples: list[str]) -> VCFGeno:
		rss = VCFFillable.collect_records(vcfs)
		sample_table = [ vcf.samples for vcf in vcfs ]
		return VCFFillable.integrate(vcfs[0].vcf, rss,
										sample_table, orig_samples)


__all__ = ['VCFFillable']
