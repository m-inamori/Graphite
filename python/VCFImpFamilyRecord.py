from __future__ import annotations

# coding: utf-8
# VCFImpFamily.py
# VCFFamilyを継承するこのプロジェクト用のクラス

from abc import ABCMeta, abstractmethod
from functools import reduce
from collections import defaultdict
from operator import add, xor
from enum import Enum

from VCFFamily import *
from TypeDeterminer import ParentComb
from common import is_all_same


#################### FillType ####################

class FillType(Enum):
	MAT = 0
	PAT = 1
	FILLED = 2
	IMPUTABLE = 3
	UNABLE = 4


#################### VCFImpFamilyRecord ####################

class VCFImpFamilyRecord(VCFFamilyRecord, metaclass=ABCMeta):
	def __init__(self, pos: int, geno: list[int], index: int,
						parents_wrong_type: str, pair: ParentComb) -> None:
		super().__init__(pos, geno)
		self.index: int = index
		self.parents_wrong_type: str = parents_wrong_type
		self.pair: ParentComb = pair
	
	def is_fixed(self) -> bool:
		return (self.pair in (ParentComb.P00x00, ParentComb.P11x11) or
				(self.pair.is_heterohomo() and self.is_right()))
	
	def is_right(self) -> bool:
		return self.parents_wrong_type == 'Right'
	
	def enable_modification(self) -> None:
		self.parents_wrong_type = 'Modifiable'
	
	# gt: 0 => 0|0 x 1|1 gt: 2 => 1|1 x 0|0
	def set_00x11_parents(self, i: int, gt: int) -> None:
		j = 1 if i == 0 else 0
		self.geno[i] = gt // 2 * 3 + 4
		self.geno[j] = 7 - gt // 2 * 3
		prog_gt = 5 if xor(gt == 0, i == 0) else 6
		for k in range(2, len(self.geno)):
			self.geno[k] = prog_gt
	
	@abstractmethod
	def is_imputable(self) -> bool:
		pass
	
	@abstractmethod
	def get_fill_type(self) -> FillType:
		pass
	
	# 複数の家系で同じサンプルを集めたRecord
	# あるGenotypeはfixedのみ、その他はfixed以外
	# -> (right GT, [(wrong record, parent index)])
	@staticmethod
	def which_is_fixed(v: list[tuple[VCFImpFamilyRecord, int]]
						) -> tuple[int, list[tuple[VCFImpFamilyRecord, int]]]:
		gts = [ r.unphased(i) for r, i in v ]
		if is_all_same(gts):
			return (gts[0], [])
		
		# collect reads in Genotype
		# ここで./.を無視してもよいはず
		dic = defaultdict(list)	# { int gt: [(record, parent index)] }
		for (r, i), gt in zip(v, gts):
			dic[gt].append((r, i))
		items = list(dic.items())	# [(int gt, [(record, parent index)])]
		
		bs = [ all(r.is_fixed() for r, i in w) for gt, w in items ]
		# 全部fixedが一つでないとあきらめる
		# 先頭のGenotypeを返すが、空のレコードのリストを返すので何でもよい
		if sum(1 for b in bs if b) != 1:
			return (gts[0], [])
		
		bs2 = [ all(not r.is_fixed() for r, i in w) for gt, w in items ]
		if sum(1 for b in bs2 if not b) != 1:
			return (gts[0], [])
		
		fixed_index = next(i for i, b in enumerate(bs) if b)
		fixed_GT = items[fixed_index][0]
		if fixed_GT == 1:	# 0/1
			return (gts[0], [])
		
		# 0/0 x 1/1パターン以外のRecordはあとで直す機会がある
		wrongs = [ (r, k) for j, item in enumerate(items) if j != fixed_index
							for r, k in item[1] if r.pair == ParentComb.P00x11 ]
		return (fixed_GT, wrongs)
	
	# 子どものGenotypeから0/0 x 1/1と推定される場合、
	# 親のどちらが0/0でどちらが1/1なのかを間違えている場合が多い
	# 親が複数の家系にまたがっているとき、修正できるかもしれない
	@staticmethod
	def modify_00x11(records: list[tuple[VCFImpFamilyRecord, int]]) -> None:
		right_gt, wrongs = VCFImpFamilyRecord.which_is_fixed(records)
		for record, i in wrongs:
			record.set_00x11_parents(i, right_gt)
