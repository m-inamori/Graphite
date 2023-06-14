from __future__ import annotations

# coding: utf-8
# VCFIsolated.py
# 孤立したサンプル

from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator

from VCFFamily import *
from VCFFillable import *
from Map import *
import Imputer
from option import *

Haplotype = Tuple[Tuple[int, int], Tuple[int, int]]


#################### VCFIsolated ####################

# imputed samples at the beginning, followed by reference samples
class VCFIsolated(VCFSmall, VCFMeasurable):
	def __init__(self, header: list[list[str]], num_imputed_samples: int,
										records: list[VCFRecord], map_: Map):
		VCFSmall.__init__(self, header, records)
		VCFMeasurable.__init__(self, map_)
		self.num_imputed_samples = num_imputed_samples
	
	def is_block(self, record: VCFRecord, rs: list[VCFRecord]) -> bool:
		length = self.cM(record.pos()) - self.cM(rs[0].pos())
		if length < 1.0:
			return True
		elif len(rs) < 10 and length < 10.0:
			return True
		else:
			return False
	
	# divide VCF by 1 cM
	def divide_by_cM(self) -> Iterator[VCFIsolated]:
		# 1cMを超えても10個以内で10cM以内なら塊とみなす
		rs = self.records[:1]
		for i in range(1, len(self)):
			record = self.records[i]
			if self.is_block(record, rs):
				rs.append(record)
			else:
				yield VCFIsolated(self.header, self.num_imputed_samples,
															rs, self.map)
				rs = [record]
		
		yield VCFIsolated(self.header, self.num_imputed_samples,
															rs, self.map)
	
	def score_each(self, hap: Haplotype, i: int, j: int) -> int:
		record = self.records[j]
		if record.is_NA(i):
			return 0
		
		hap1, hap2 = hap
		i1, j1 = hap1
		i2, j2 = hap2
		int_gt = record.get_int_gt(i)
		h1 = int(record.v[i1+9][j1*2])
		h2 = int(record.v[i2+9][j2*2])
		return 1 if int_gt == h1 + h2 else 0
	
	def score(self, hap: Haplotype, i: int) -> int:
		total_score = 0
		for j in range(len(self.records)):
			total_score += self.score_each(hap, i, j)
		return total_score
	
	def generate_haplotype_pairs(self):
		range_obj = range(self.num_imputed_samples, len(self.samples))
		for i1, i2 in product(range_obj, repeat=2):
			if i1 == i2:
				continue
			for j1, j2 in product(range(2), repeat=2):	# which haplotype
				yield (i1, j1), (i2, j2)
	
	def collect_optimal_haplotype_pairs(self, i: int) -> list[Haplotype]:
		max_score = 0
		max_combs: list[Haplotype] = []
		for hap in self.generate_haplotype_pairs():
			s = self.score(hap, i)
			if s > max_score:
				max_score = s
				max_combs = [hap]
			elif s == max_score:
				max_combs.append(hap)
		return max_combs
	
	def set_haplotype(self, hap: Haplotype, i: int):
		hap1, hap2 = hap
		i1, j1 = hap1
		i2, j2 = hap2
		for record in self.records:
			gt1 = record.v[i1+9][j1*2]
			gt2 = record.v[i2+9][j2*2]
			record.set_GT(i, gt1 + '|' + gt2)
	
	def impute_cM_each_sample(self, prev_hap: Haplotype, i: int) -> Haplotype:
		def score(hap: Haplotype, prev_hap: Haplotype) -> int:
			hap1, hap2 = hap
			prev_hap1, prev_hap2 = prev_hap
			return ((1 if hap1 == prev_hap1 else 0) +
					(1 if hap2 == prev_hap2 else 0))
		
		# とりあえず、総当たりにしてみる
		combs = self.collect_optimal_haplotype_pairs(i)
		
		# 前との一致度を計算する
		v = [ (score(hap, prev_hap), hap) for hap in combs ]
		max_score, _ = max(v)
		filtered_combs = [ hap for s, hap in v if s == max_score ]
		# 複数の候補があれば乱数っぽく決める
		j = self.records[0].pos() % len(filtered_combs)
		hap = filtered_combs[j]
		self.set_haplotype(hap, i)
		return hap
	
	def impute_cM(self, prev_haps: list[Haplotype]) -> list[Haplotype]:
		haps: list[Haplotype] = []
		for i, prev_hap in enumerate(prev_haps):
			hap = self.impute_cM_each_sample(prev_hap, i)
			haps.append(hap)
		return haps
	
	def impute(self):
		haps = [((0, -1), (0, -1))] * self.num_imputed_samples
		for vcf_cM in self.divide_by_cM():
			haps = vcf_cM.impute_cM(haps)
	
	def extract_isolated_samples(self) -> VCFSmall:
		isolated_samples = self.samples[:self.num_imputed_samples]
		return self.extract_samples(isolated_samples)
	
	@staticmethod
	def create(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
						samples: list[str], references: list[str],
						gmap: Map, num_threads: int) -> list[VCFIsolated]:
		sample_columns = orig_vcf.extract_columns(samples)
		ref_columns = merged_vcf.extract_columns(references)
		column_table = VCFIsolated.divide_columns(sample_columns, num_threads)
		vcfs: list[VCFIsolated] = []
		for cs in column_table:
			div_samples = [ orig_vcf.samples[c-9] for c in cs ]
			new_samples = div_samples + references
			header = orig_vcf.trim_header(new_samples)
			records: list[VCFRecord] = []
			for j in range(len(orig_vcf)):
				record = orig_vcf.records[j]
				v = (record.v[:9] +
						[ record.v[c] for c in cs ] +
						[ merged_vcf.records[j].v[c] for c in ref_columns ])
				new_record = VCFRecord(v, new_samples)
				records.append(new_record)
			vcf = VCFIsolated(header, len(cs), records, gmap)
			vcfs.append(vcf)
		return vcfs
	
	@staticmethod
	def divide_columns(cs: list[int], num: int) -> list[list[int]]:
		if len(cs) <= num:
			return [ [c] for c in cs ]
		
		css: list[list[int]] = [ [] for _ in range(num) ]
		for i, c in enumerate(cs):
			css[i%num].append(c)
		return css

__all__ = ['VCFIsolated']
