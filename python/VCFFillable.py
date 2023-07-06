from __future__ import annotations

# coding: utf-8

from functools import reduce
from itertools import *
from collections import Counter, defaultdict
from math import log

from typing import Optional

from VCFFamily import *
from VCFImpFamily import VCFImpFamilyRecord
from VCFHeteroHomo import *
from option import *
from common import *


#################### Genotype ####################

class Genotype:
	def __init__(self, s: str):
		self.phasing: bool	= s[1] == '|'
		self.gt1: str		= s[0]
		self.gt2: str		= s[2]
	
	def gts(self) -> list[str]:
		return [self.gt1, self.gt2]
	
	def __str__(self) -> str:
		return self.gt1 + ('|' if self.phasing else '/') + self.gt2
	
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
	def is_valid(gt: str, mat_gt: int, pat_gt: int) -> bool:
		mat_gts = Genotype.possible_gts(mat_gt)
		pat_gts = Genotype.possible_gts(pat_gt)
		return len(gt) >= 3 and (gt[0] in mat_gts and gt[2] in pat_gts or
								 gt[2] in mat_gts and gt[0] in pat_gts)
	
	# ここにあるのがいいのかわからない
	# できれば違うファイルにして、これ関連をここに集めたい
	@staticmethod
	def possible_gts(gt: int) -> list[str]:
		if gt == 0:
			return ['0']
		elif gt == 3:
			return ['1']
		else:
			return ['0', '1']


#################### VCFFillableRecord ####################

class VCFFillableRecord(VCFFamilyRecord):
	def __init__(self, v: list[str], samples: list[str], index: int,
														type: str, pair: int):
		super().__init__(v, samples)
		self.index: int = index
		self.type: str = type
		self.pair: int = pair
	
	def possible_phasings(self) -> list[tuple[int, int]]:
		# pair 0: 0x0, 1: 0x1, 2: 0x2, 3: 1x1, 4: 1x2, 5: 2x2
		if self.pair == 0:		# 0/0 x 0/0
			return [(0, 0)]
		elif self.pair == 1:	# 0/0 x 0/1
			return [(0, 1), (0, 2), (1, 0), (2, 0),
					(1, 1), (1, 2), (2, 1), (2, 2)]
		elif self.pair == 2:	# 0/1 x 0/1
			return [(1, 1), (1, 2), (2, 1), (2, 2),
					(0, 1), (0, 2), (1, 0), (2, 0),
					(1, 3), (2, 3), (3, 1), (3, 2)]
		elif self.pair == 3:	# 0/0 x 1/1
			return [(0, 3), (3, 0)]
		elif self.pair == 4:	# 0/1 x 1/1
			return [(1, 3), (2, 3), (3, 1), (3, 2),
					(1, 1), (1, 2), (2, 1), (2, 2)]
		elif self.pair == 5:	# 1/1 x 1/1
			return [(3, 3)]
		else:
			# all
			return list(product(range(4), range(4)))
	
	def group_ids(self) -> list[int]:
		w = self.v[7].split('=')
		return list(map(int, w[1].split(',')))
	
	def set_group_ids(self, ids: list[int]):
		self.v[7] = 'GR=%s' % (','.join(map(str, ids)))
	
	def gt_from_parent(self, mat_from: int, pat_from: int) -> str:
		gt_from_mat = '.' if mat_from == 0 else self.v[9][mat_from*2-2]
		gt_from_pat = '.' if pat_from == 0 else self.v[10][pat_from*2-2]
		return gt_from_mat + '|' + gt_from_pat
	
	def gt_from_mat(self, mat_from: int, c: int) -> str:
		gt = self.get_GT(c - 9)
		int_gt = self.get_int_gt(c - 9)
		mat_gt = self.v[9][mat_from*2-2]
		pat_GT = self.pat_GT()
		if int_gt == 0:
			if mat_gt == '0':
				return '0|0'
			elif '0' in pat_GT:
				return '1|0'
			else:
				return '1|1'
		elif int_gt == 2:
			if mat_gt == '1':
				return '1|1'
			elif '1' in pat_GT:
				return '0|1'
			else:
				return '0|0'
		elif int_gt == 1:
			if mat_gt == '0':
				return '0|1'
			else:
				return '1|0'
		else:
			return '.|.'
	
	def gt_from_pat(self, pat_from: int, c: int) -> str:
		gt = self.get_GT(c - 9)
		int_gt = self.get_int_gt(c - 9)
		pat_gt = self.v[10][pat_from*2-2]
		mat_GT = self.mat_GT()
		if int_gt == 0:
			if pat_gt == '0':
				return '0|0'
			elif '0' in mat_GT:
				return '0|1'
			else:
				return '1|1'
		elif int_gt == 2:
			if pat_gt == '1':
				return '1|1'
			elif '1' in mat_GT:
				return '1|0'
			else:
				return '0|0'
		elif int_gt == 1:
			if pat_gt == '0':
				return '1|0'
			else:
				return '0|1'
		else:
			return '.|.'
	
	def is_mat_hetero(self) -> bool:
		return self.v[9][0] != self.v[9][2]
	
	def mat_from(self, c: int) -> int:
		if self.v[c][0] == '.':
			return 0
		elif not self.is_hetero(0):
			return 0
		elif self.v[c][0] == self.v[9][0]:
			return 1
		else:
			return 2
	
	def pat_from(self, c: int) -> int:
		if self.v[c][0] == '.':
			return 0
		elif not self.is_hetero(1):
			return 0
		elif self.v[c][2] == self.v[10][0]:
			return 1
		else:
			return 2
	
	def find_geno_type(self, type: str) -> int:
		for i, t in enumerate(self.v[8].split(':')):
			if t == type:
				return i
		else:
			return -1
	
	def fill_PGT(self):
		i_GT: int = self.find_geno_type('GT')
		assert(i_GT != -1)
		i_PGT: int = self.find_geno_type('PGT')
		if i_PGT == -1:
			return
		
		for j in range(9, len(self.v)):
			v = self.v[j].split(':')
			if i_PGT >= len(v):		# 補完の際にGTだけになった
				continue
			v[i_PGT] = v[i_GT]
			self.v[j] = ':'.join(v)
	
	def modify_parents_type(self):
		if (self.pair != 3 and
				((self.v[9][:3] == '0|0' and self.v[10][:3] == '1|1') or
				 (self.v[9][:3] == '1|1' and self.v[10][:3] == '0|0'))):
			self.pair = 3
	
	@staticmethod
	def convert(record: VCFImpFamilyRecord) -> VCFFillableRecord:
		type = record.get_fill_type()
		return VCFFillableRecord(record.v, record.samples, record.index,
															type, record.pair)
	
	@staticmethod
	def merge(records: list[VCFFillableRecord],
								samples: list[str]) -> VCFRecord:
		v = records[0].v[:]
		for record in records[1:]:
			v.extend(record.v[9:])
		return VCFRecord(v, samples)
	
	def decide_by_majority(self, GTs: list[str]) -> str:
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
			i = reduce(lambda x, y: x + int(y), self.v[1], 0) % d
			return candidate_GTs[i]
	
	def swap_parents(self, pos: int, GT: str):
		if GT == self.v[pos+9][:3]:
			return
		elif GT in ('0|0', '1|1'):
			is_mat_00 = (pos == 0) ^ (GT == '1|1')
			self.set_GT(0, '0|0' if is_mat_00 else '1|1')
			self.set_GT(1, '1|1' if is_mat_00 else '0|0')
			# 子どもも入れ換えなければならない
			prog_GT = '0|1' if is_mat_00 else '1|0'
			for i in range(2, len(self.samples)):
				self.set_GT(i, prog_GT)
	
	@staticmethod
	def decide_duplicated_Genotype(records: list[VCFFillableRecord],
											ps: list[tuple[int, int]]) -> str:
		# 全部./.ならそのまま
		# 子どもがいたらそれ
		# ./.を除いて全部同じだったらそれ
		# ./.と0/0 x 1/1で全部なら0/0 x 1/1で多数決
		# ./.と0/0 x 1/1を除いて全部同じだったらそれ
		# ./.と0/0 x 1/1を除いて複数種あれば多数決
		GTs = [ records[i].v[j+9][:3] for i, j in ps ]
		if all(GT == './.' for GT in GTs):
			return './.'
		
		# 子どもが優先
		for (i, j), GT in zip(ps, GTs):
			if j >= 2 and GT != './.':
				return GT
		
		GTs2 = [ GT for GT in GTs if GT != './.' ]
		if is_all_same(GTs2):
			return GTs2[0]
		
		GTs3 = [ GT for (i, j), GT in zip(ps, GTs)
								if GT != './.' and records[i].pair != 3 ]
		if not GTs3:	# 0/0 x 1/1しかない
			return records[0].decide_by_majority(GTs2)
		elif is_all_same(GTs3):
			return GTs3[0]
		else:
			return records[0].decide_by_majority(GTs3)
	
	@staticmethod
	def integrate_each_sample(records: list[VCFFillableRecord],
												ps: list[tuple[int, int]]):
		GT = VCFFillableRecord.decide_duplicated_Genotype(records, ps)
		
		for i, j in ps:
			# 0/0 x 1/1の親なら相手と入れ替えなければならない
			if records[i].pair == 3 and j <= 1:
				records[i].swap_parents(j, GT)
	
	@staticmethod
	def integrate(records: list[VCFFillableRecord], samples: list[str],
						pos_samples: list[list[tuple[int, int]]]) -> VCFRecord:
		# 0|0 x 1|1ならひっくり返せる
		for ps in pos_samples:
			GTs = [ records[i].v[j+9][:3] for i, j in ps ]
			if not is_all_same(GTs):
				VCFFillableRecord.integrate_each_sample(records, ps)
		
		# ひっくり返した後にGTをまとめる
		v: list[str] = records[0].v[:9]
		for ps in pos_samples:
			i, j = ps[0]
			v.append(records[i].v[j+9])
		return VCFRecord(v, samples)


#################### RecordSet ####################

class RecordSet:
	def __init__(self, r, pm, nm, pp, np):
		self.record: Optional[VCFFillableRecord] = r
		self.prev_mat_record: Optional[VCFFillableRecord] = pm
		self.next_mat_record: Optional[VCFFillableRecord] = nm
		self.prev_pat_record: Optional[VCFFillableRecord] = pp
		self.next_pat_record: Optional[VCFFillableRecord] = np
	
	def records(self) -> list[Optional[VCFFillableRecord]]:
		return [self.record, self.prev_mat_record, self.next_mat_record,
							 self.prev_pat_record, self.next_pat_record]
	
	def gt_each(self, i: int, record: Optional[VCFFillableRecord]) -> str:
		return "" if record is None else record.v[i+9]
	
	def gt(self, i: int) -> str:
		return self.gt_each(i, self.record)
	
	def __select_phasing(self, candidates: list[tuple[int, int]]
												) -> tuple[int, int]:
		if len(candidates) == 1 or self.record is None:
			return candidates[0]
		
		mat_int_gt = self.record.mat_int_gt()
		pat_int_gt = self.record.pat_int_gt()
		
		v: list[tuple[float, int, int]] = []
		for mat_phasing, pat_phasing in candidates:
			mat_int_gt1 = (mat_phasing >> 1) + (mat_phasing & 1)
			pat_int_gt1 = (pat_phasing >> 1) + (pat_phasing & 1)
			score = (abs(mat_int_gt1 - mat_int_gt) +
					 abs(pat_int_gt1 - pat_int_gt))
			v.append((score, mat_phasing, pat_phasing))
		
		# これだと、同じスコアならphasingが小さい方から選んでいる
		# 本当はランダム的に選びたい
		_, mat_phasing, pat_phasing = min(v)
		return (mat_phasing, pat_phasing)
	
	def determine_phasing_core(self, lls: list[tuple[float, int, int]]
														) -> tuple[int, int]:
		# 最大に近いllを集める
		candidates = []
		max_ll = lls[-1][0]
		for ll, mat_phasing, pat_phasing in reversed(lls):
			if ll > max_ll - 1e-9:
				candidates.append((mat_phasing, pat_phasing))
			else:
				break
		return self.__select_phasing(candidates)


#################### VCFFillable ####################

class VCFFillable(VCFBase, VCFSmallBase, VCFFamilyBase):
	def __init__(self, header: list[list[str]],
							records: list[VCFFillableRecord]):
		VCFBase.__init__(self, header)
		VCFFamilyBase.__init__(self)
		self.records: list[VCFFillableRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def modify(self):
		# typeが'UNABLE', 'IMPUTABLE', 'MAT', 'PAT', 'FIXED'でrecordを分ける
		groups: list[tuple[str, list[VCFFillableRecord]]] = [ (g, list(v))
					for g, v in groupby(self.records, key=lambda r: r.type) ]
		for i, (key, subrecords) in enumerate(groups):
			if key == 'IMPUTABLE' or key == 'UNABLE':
				self.__phase(groups, i, True)
		
		for i, record in enumerate(self.records):
			# 家系ごとで./.にしたGenotypeを補完
			if record.type == 'MAT':
				self.__impute_NA_mat(i)
			elif record.type == 'PAT':
				self.__impute_NA_pat(i)
			elif record.type in ('IMPUTABLE', 'UNABLE'):
				self.__impute_others(i)
		
		for record in self.records:
			record.fill_PGT()
	
	def phase_hetero_hetero(self):
		# typeが'IMPUTABLE', 'MAT', 'PAT', 'FIXED'でrecordを分ける
		groups: list[tuple[str, list[VCFFillableRecord]]] = [ (g, list(v))
					for g, v in groupby(self.records, key=lambda r: r.type) ]
		for i, (key, subrecords) in enumerate(groups):
			if key == 'IMPUTABLE':
				self.__phase(groups, i, False)
		
		# この部分はフルで必要？
		for i, record in enumerate(self.records):
			# 家系ごとで./.にしたGenotypeを補完
			if record.type == 'MAT':
				self.__impute_NA_mat(i)
			elif record.type == 'PAT':
				self.__impute_NA_pat(i)
			elif record.type in ('IMPUTABLE', 'UNABLE'):
				self.__impute_others(i)
	
	def __phase(self, groups: list[tuple[str, list[VCFFillableRecord]]],
											i: int, necessary_parents_phasing):
		prev_mat_record = self.find_prev_record(groups, i, 'MAT')
		next_mat_record = self.find_next_record(groups, i, 'MAT')
		prev_pat_record = self.find_prev_record(groups, i, 'PAT')
		next_pat_record = self.find_next_record(groups, i, 'PAT')
		key, records = groups[i]
		for record in records:
			record_set = RecordSet(record, prev_mat_record,
							next_mat_record, prev_pat_record, next_pat_record)
			if necessary_parents_phasing:
				self.__determine_parents_phasing(record_set)
			self.__impute_core(record_set)
	
	def find_prev_record(self,
							groups: list[tuple[str, list[VCFFillableRecord]]],
							i: int, g: str) -> Optional[VCFFillableRecord]:
		key, records = groups[i]
		chr = records[0].v[0]
		for j in range(i-1, -1, -1):
			key, records = groups[j]
			if records[0].v[0] != chr:
				return None
			if key == g:
				return records[-1]
		else:
			return None
	
	def find_next_record(self,
							groups: list[tuple[str, list[VCFFillableRecord]]],
							i: int, g: str) -> Optional[VCFFillableRecord]:
		key, records = groups[i]
		chr = records[0].v[0]
		for j in range(i+1, len(groups)):
			key, records = groups[j]
			if records[0].v[0] != chr:
				return None
			if key == g:
				return records[0]
		else:
			return None
	
	def __determine_parents_phasing(self, recordset: RecordSet):
		if recordset.record is None:
			return
		
		record = recordset.record
		lls = []
		for mat_p, pat_p in record.possible_phasings():
			ll = self.compute_phasing_likelihood(recordset, mat_p, pat_p)
			lls.append((ll, mat_p, pat_p))
		lls.sort()
		
		mat_gt, pat_gt = recordset.determine_phasing_core(lls)
		if record is None:
			return
		gt = ['0|0', '1|0', '0|1', '1|1']
		record.v[9] = gt[mat_gt] + record.v[9][3:]
		record.v[10] = gt[pat_gt] + record.v[10][3:]
	
	# mat_gt, pat_gt : 0|0 0|1 1|0 1|1を0～3で表す
	def compute_phasing_likelihood(self, recordset: RecordSet,
								mat_gt: int, pat_gt: int) -> float:
		def gts(record: Optional[VCFFillableRecord]) -> list[str]:
			return record.v[11:] if record else [''] * (len(self.samples) - 2)
		
		memo = { (0, 0): [0.5, 0.5], (0, 1): [0.9, 0.1], (0, 2): [0.1, 0.9],
				 (1, 0): [0.9, 0.1], (1, 1): [0.99, 0.01], (1, 2): [0.5, 0.5],
				 (2, 0): [0.1, 0.9], (2, 1): [0.5, 0.5], (2, 2): [0.01, 0.99] }
		def probs_from_which_chrom(prev_chrom: int,
								   next_chrom: int) -> list[float]:
			return memo[(prev_chrom, next_chrom)]
		
		def likelihood_each(gt: str, probs_mat: list[float],
									 probs_pat: list[float]) -> float:
			sum_gt = int(gt[0]) + int(gt[2])
			likelihood = sum(probs_mat[i] * probs_pat[j]
								for i, j in product(range(2), repeat=2)
								if ((mat_gt >> i) & 1) +
								   ((pat_gt >> j) & 1) == sum_gt)
			return log(likelihood)
		
		record, mat_record1, mat_record2, pat_record1, pat_record2 = \
															recordset.records()
		ll = 0.0	# log of likelihood
		for gt, mat_gt1, mat_gt2, pat_gt1, pat_gt2 in zip(*map(gts,
														recordset.records())):
			if not Genotype.is_valid(gt, mat_gt, pat_gt):
				ll += log(0.0001)
				continue
			prev_mat_from = self.from_which_chrom(mat_gt1, mat_record1, True)
			next_mat_from = self.from_which_chrom(mat_gt2, mat_record2, True)
			prev_pat_from = self.from_which_chrom(pat_gt1, pat_record1, False)
			next_pat_from = self.from_which_chrom(pat_gt2, pat_record2, False)
			probs_mat = probs_from_which_chrom(prev_mat_from, next_mat_from)
			probs_pat = probs_from_which_chrom(prev_pat_from, next_pat_from)
			ll += likelihood_each(gt, probs_mat, probs_pat)
		return ll
	
	def __impute_core(self, recordset: RecordSet):
		def gts(record: Optional[VCFFillableRecord]) -> list[str]:
			return record.v[11:] if record else [''] * (len(self.samples) - 2)
		
		def sum_gt(gt: str) -> int:
			try:
				return int(gt[0]) + int(gt[2])
			except ValueError:
				return -1
		
		def both_nearests() -> tuple[int, int]:
			# ここに来るときは、全てNoneではない
			if (record is None or mat_record1 is None or mat_record2 is None
							   or pat_record1 is None or pat_record2 is None):
				return (0, 0)
			if record.pos() * 2 < mat_record1.pos() + mat_record2.pos():
				mat_from = prev_mat_from
			else:
				mat_from = next_mat_from
			if record.pos() * 2 < pat_record1.pos() + pat_record2.pos():
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
				if record.pos() * 2 < pat_record1.pos() + pat_record2.pos():
					return (pairs[0][0], prev_pat_from)
				else:
					return (pairs[0][0], next_pat_from)
			elif pairs[0][1] == pairs[1][1]:	# patが同じ
				if record is None or mat_record1 is None or mat_record2 is None:
					return (0, 0)
				if record.pos() * 2 < mat_record1.pos() + mat_record2.pos():
					return (prev_mat_from, pairs[0][1])
				else:
					return (next_mat_from, pairs[0][1])
			else:	# 両親とも乗り換えている（滅多にない）
				return both_nearests()
		
		def select_pair(pairs: list[tuple[int, int]], gt: str,
									selected: bool = False) -> tuple[int, int]:
			if record is None:
				return (0, 0)
			elif not pairs:
				return (0, 0)
			elif len(pairs) == 1:
				return pairs[0]
			elif gt[0] == '.':
				return nearest_froms(pairs)
			elif selected:
				return nearest_froms(pairs)
			else:
				new_pairs = [ v for v in pairs
							if sum_gt(gt) == sum_gt(record.gt_from_parent(*v)) ]
				pair = select_pair(new_pairs, gt, True)
				if pair != (0, 0):
					return pair
				else:
					return select_pair(pairs, gt, True)
		
		def is_same_gts(gt1: str, gt2: str) -> bool:
			if gt2 == '0/1':
				return gt1 == '0|1' or gt1 == '1|0'
			elif gt2 == '0/0':
				return gt1 == '0|0'
			elif gt2 == '1/1':
				return gt1 == '1|1'
			else:
				return False
		
		def is_near(gts1: list[str], gts2: list[str]) -> bool:
			num = sum(1 for gt2 in gts2 if gt2 != '0/1')
			dist = sum(1 for gt1, gt2 in zip(gts1, gts2)
										if not is_same_gts(gt1, gt2))
			return dist < max(1, num // 2)
		
		def inverse_gts(gts: list[str], inv_mat: bool,
										inv_pat: bool) -> list[str]:
			def inverse(s: str) -> str:
				return '1' if s == '0' else '0'
			
			if inv_mat:
				if inv_pat:
					return [ inverse(gt[0]) + '|' + inverse(gt[2]) + gt[3:]
																for gt in gts ]
				else:
					return [ inverse(gt[0]) + '|' +  gt[2:] for gt in gts ]
			else:
				if inv_pat:
					return [ gt[:1] + '|' + inverse(gt[2]) + gt[3:] for gt in gts ]
				else:
					return [ gt[:1] + '|' + gt[2:] for gt in gts ]
		
		def inverse_parents_gts(record, inv_mat, inv_pat):
			# both must be hetero
			if inv_mat:
				gt = record.v[9]
				record.v[9] = gt[2] + '|' + gt[0] + gt[3:]
			if inv_pat:
				gt = record.v[10]
				record.v[10] = gt[2] + '|' + gt[0] + gt[3:]
		
		def modify_gts(gts: list[str], record: Optional[VCFFillableRecord]):
			if record is None:
				return
			
			orig_gts = record.v[11:]
			for inv_mat, inv_pat in product((False, True), repeat=2):
				new_gts = inverse_gts(gts, inv_mat, inv_pat)
				if is_near(new_gts, orig_gts):
					inverse_parents_gts(record, inv_mat, inv_pat)
					for c in range(11, len(record.v)):
						record.v[c] = new_gts[c-11]
					break
			else:
				new_gts = gts
			for c in range(11, len(record.v)):
				record.v[c] = new_gts[c-11]
		
		# どちらから来たか決める
		def select_from(froms: list[int], record1: Optional[VCFRecord],
										  record2: Optional[VCFRecord]) -> int:
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
				elif record.pos() * 2 < record1.pos() + record2.pos():
					return froms[0]
				else:
					return froms[1]
		
		def modify_parents_type(record):
			if record is not None:
				record.modify_parents_type()
		
		record, mat_record1, mat_record2, pat_record1, pat_record2 = \
															recordset.records()
		if record is None:
			return
		new_gts = record.v[11:]
		for c, gt, mat_gt1, mat_gt2, pat_gt1, pat_gt2 in zip(count(11),
												*map(gts, recordset.records())):
			prev_mat_from = self.from_which_chrom(mat_gt1, mat_record1, True)
			next_mat_from = self.from_which_chrom(mat_gt2, mat_record2, True)
			prev_pat_from = self.from_which_chrom(pat_gt1, pat_record1, False)
			next_pat_from = self.from_which_chrom(pat_gt2, pat_record2, False)
			if ((prev_mat_from == 0 and next_mat_from == 0) or
					(prev_pat_from == 0 and next_pat_from == 0)):
				new_gts[c-11] = record.v[c]
				continue
			
			mat_froms = unique_list(prev_mat_from, next_mat_from)
			pat_froms = unique_list(prev_pat_from, next_pat_from)
			pairs = [ x for x in product(mat_froms, pat_froms)
										if x[0] != 0 and x[1] != 0 ]
			mat_from, pat_from = select_pair(pairs, gt)
			res = record.v[c][3:]
			if (mat_from, pat_from) == (0, 0):
				if any(mat_from != 0 for mat_from in mat_froms):
					mat_from_nz = select_from(mat_froms,
												mat_record1, mat_record2)
					new_gts[c-11] = record.gt_from_mat(mat_from_nz, c) + res
				elif any(pat_from != 0 for pat_from in pat_froms):
					pat_from_nz = select_from(pat_froms,
												pat_record1, pat_record2)
					new_gts[c-11] = record.gt_from_pat(pat_from_nz, c) + res
				else:
					continue	# phasingしない
			else:
				new_gts[c-11] = record.gt_from_parent(mat_from, pat_from) + res
		
		modify_gts(new_gts, record)		# 両親が間違っていた時の処理
		modify_parents_type(record)
	
	# phasingされている前提
	def from_which_chrom(self, gt: str, record: Optional[VCFFillableRecord],
															mat: bool) -> int:
		if record is None:
			return 0
		
		i = 0 if mat else 1
		parent_gt = record.v[9+i]
		return 1 if parent_gt[0] == gt[i*2] else 2
	
	def find_prev_same_type_record(self, i: int, c: int
											) -> Optional[VCFFillableRecord]:
		type = self.records[i].type
		chromosome = self.records[i].v[0]
		for j in range(i - 1, -1, -1):
			record = self.records[j]
			if record.v[0] != chromosome:
				return None
			if record.type == type and record.v[c][:3] != './.':
				return self.records[j]
		else:
			return None
	
	def find_next_same_type_record(self, i: int, c: int
											) -> Optional[VCFFillableRecord]:
		type = self.records[i].type
		chromosome = self.records[i].v[0]
		for j in range(i + 1, len(self.records)):
			record = self.records[j]
			if record.v[0] != chromosome:
				return None
			if record.type == type and record.v[c][:3] != './.':
				return self.records[j]
		else:
			return None
	
	# 家系ごとで./.にしたGenotypeを補完
	def __impute_NA_mat(self, i: int):
		record = self.records[i]
		for c in range(11, len(record.v)):
			if record.v[c][:3] == '.|.':
				self.__impute_NA_mat_each(i, c)
	
	def __impute_NA_mat_each(self, i: int, c: int):
		def select_mat(pairs: list[tuple[int, int]]) -> int:
			if len(pairs) == 1:
				return pairs[0][0]
#			elif prev_record is None or next_record is None:
			elif prev_record is None and next_record is None:
				return 0
			elif next_record is None:
				return prev_mat_from
			elif prev_record is None:
				return next_mat_from
			elif record.pos() * 2 < prev_record.pos() + next_record.pos():
				return prev_mat_from
			else:
				return next_mat_from
		
		record = self.records[i]
		prev_record = self.find_prev_same_type_record(i, c)
		next_record = self.find_next_same_type_record(i, c)
		prev_gt = prev_record.v[c] if prev_record is not None else ''
		next_gt = next_record.v[c] if next_record is not None else ''
		prev_mat_from = self.from_which_chrom(prev_gt, prev_record, True)
		next_mat_from = self.from_which_chrom(next_gt, next_record, True)
		mat_froms = unique_list(prev_mat_from, next_mat_from)
		pairs = [ (x, 1) for x in mat_froms if x != 0 ]
		if not pairs:	# 両側ともにmatが無い（まず無い）
			return
		
		mat_from = select_mat(pairs)
		if mat_from == 0:
			return
		
		record.v[c] = record.gt_from_parent(mat_from, 1) + record.v[c][3:]
	
	# 家系ごとで./.にしたGenotypeを補完
	def __impute_NA_pat(self, i: int):
		record = self.records[i]
		for c in range(11, len(record.v)):
			if record.v[c][:3] == './.':
				self.__impute_NA_pat_each(i, c)
	
	def __impute_NA_pat_each(self, i: int, c: int):
		def select_pat(pairs: list[tuple[int, int]]) -> int:
			if len(pairs) == 1:
				return pairs[0][1]
#			elif prev_record is None or next_record is None:
			elif prev_record is None and next_record is None:
				return 0
			elif next_record is None:
				return prev_pat_from
			elif prev_record is None:
				return next_pat_from
			elif record.pos() * 2 < prev_record.pos() + next_record.pos():
				return prev_pat_from
			else:
				return next_pat_from
		
		record = self.records[i]
		prev_record = self.find_prev_same_type_record(i, c)
		next_record = self.find_next_same_type_record(i, c)
		prev_gt = prev_record.v[c] if prev_record is not None else ''
		next_gt = next_record.v[c] if next_record is not None else ''
		prev_pat_from = self.from_which_chrom(prev_gt, prev_record, False)
		next_pat_from = self.from_which_chrom(next_gt, next_record, False)
		pat_froms = unique_list(prev_pat_from, next_pat_from)
		pairs = [ (1, x) for x in pat_froms if x != 0 ]
		if not pairs:	# 両側ともにpatが無い（まず無い）
			return
		
		pat_from = select_pat(pairs)
		record.v[c] = record.gt_from_parent(1, pat_from) + record.v[c][3:]
	
	def __impute_others(self, i: int):
		def correct(GT: str, c: int, record: VCFFillableRecord) -> str:
			int_gt = record.get_int_gt(c-9)
			if int_gt == -1:
				if GT[0] != '.':
					return GT[:2] + str(record.pos() % 2)
				elif GT[2] != '.':
					return str(record.pos() % 2) + GT[1:]
				else:
					return (str(record.pos() % 2) + '|' +
									str((record.pos() >> 1) % 2))
			elif int_gt == 0:
				if GT[0] == '1':
					return '1|0'
				elif GT[2] == '1':
					return '0|1'
				else:
					return '0|0'
			elif int_gt == 1:
				if GT[0] == '1':
					return '1|0'
				elif GT[2] == '1':
					return '0|1'
				else:
					# どちらにすべきか判別がつかないので適当に選ぶ
					return '0|1' if int(record.v[1]) % 2 == 0 else '1|0'
			else:
				if GT[0] == '0':
					return '0|1'
				elif GT[2] == '0':
					return '1|0'
				else:
					return '1|1'
		
		record = self.records[i]
		cs = [ c for c in range(11, len(record.v))
				if record.v[c][1] == '/' or record.get_int_gt(c-9) == -1 ]
		if not cs:
			return
		
		mat_homo = record.is_homo(0)
		pat_homo = record.is_homo(1)
		for c in cs:
			mat_from = 1 if mat_homo else self.__find_mat_from(i, c)
			pat_from = 1 if pat_homo else self.__find_pat_from(i, c)
			GT = record.gt_from_parent(mat_from, pat_from)
			if '.' in GT:
				GT = correct(GT, c, record)
			record.set_GT(c-9, GT)
	
	def __find_mat_from(self, i: int, c: int) -> int:
		return self.__select_from(self.__find_prev_mat_from(i, c),
									self.__find_next_mat_from(i, c), i)
	
	def __find_pat_from(self, i: int, c: int) -> int:
		return self.__select_from(self.__find_prev_pat_from(i, c),
									self.__find_next_pat_from(i, c), i)
	
	def __select_from(self, f1: tuple[int, int],
							f2: tuple[int, int], i: int) -> int:
		i1, from1 = f1
		i2, from2 = f2
		if from1 == 0 and from2 == 0:
			# 前後がないとき乱数的に決める
			r0 = self.records[i]
			return r0.pos() % 2 + 1
		
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
			if r0.pos() * 2 <= r2.pos() + r1.pos():
				return from1
			else:
				return from2
	
	def __find_prev_mat_from(self, i: int, c: int) -> tuple[int, int]:
		for k in range(i - 1, -1, -1):
			from1 = self.records[k].mat_from(c)
			if from1 != 0:
				return (k, from1)
		return (-1, 0)
	
	def __find_next_mat_from(self, i: int, c: int) -> tuple[int, int]:
		for k in range(i + 1, len(self.records)):
			from1 = self.records[k].mat_from(c)
			if from1 != 0:
				return (k, from1)
		return (-1, 0)
	
	def __find_prev_pat_from(self, i: int, c: int) -> tuple[int, int]:
		for k in range(i - 1, -1, -1):
			from1 = self.records[k].pat_from(c)
			if from1 != 0:
				return (k, from1)
		return (-1, 0)
	
	def __find_next_pat_from(self, i: int, c: int) -> tuple[int, int]:
		for k in range(i + 1, len(self.records)):
			from1 = self.records[k].pat_from(c)
			if from1 != 0:
				return (k, from1)
		return (-1, 0)
	
	def remove_parents(self) -> VCFSmall:
		samples: list[str] = self.samples[2:]
		records: list[VCFRecord] = []
		for record in self.records:
			v = record.v[:9] + record.v[11:]
			r = VCFRecord(v, samples)
			records.append(r)
		header = self.trim_header(samples)
		return VCFSmall(header, records)
	
	@staticmethod
	def fill(vcfs: list[VCFHeteroHomo],
					records: list[VCFImpFamilyRecord]) -> VCFFillable:
		merged_records = VCFFillable.merge_records(vcfs, records)
		vcf: VCFFillable = VCFFillable(vcfs[0].header, merged_records)
		vcf.modify()
		return vcf
	
	@staticmethod
	def merge_records(vcfs: list[VCFHeteroHomo],
						records: list[VCFImpFamilyRecord],
						) -> list[VCFFillableRecord]:
		all_records = [ VCFFillableRecord.convert(record) for record in
						chain((r for vcf in vcfs for r in vcf.records), records)
		]
		all_records.sort(key=lambda record: record.index)
		return all_records
	
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
	def integrate(vcf: VCFFillable, rss: list[list[VCFFillableRecord]],
										orig_samples: list[str]) -> VCFSmall:
		sss = [ r.samples for r in rss[0] ]
		samples, pos_samples = VCFFillable.integrate_samples(sss, orig_samples)
		records: list[VCFRecord] = []
		for rs in rss:
			records.append(VCFFillableRecord.integrate(rs, samples, pos_samples))
		header = vcf.trim_header(records[0].samples)
		new_vcf = VCFSmall(header, [])
		new_vcf.records.extend(records)
		return new_vcf
	
	@staticmethod
	def merge(vcfs: list[VCFFillable], orig_samples: list[str]) -> VCFSmall:
		rss = VCFFillable.collect_records(vcfs)
		return VCFFillable.integrate(vcfs[0], rss, orig_samples)


__all__ = ['VCFFillableRecord', 'VCFFillable']
