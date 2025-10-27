from __future__ import annotations

# coding: utf-8
# VCFHeteroHomo.py

from collections import defaultdict, Counter
from math import log10
from typing import List, Tuple, Optional, IO, Dict, Iterator
import random
import time

from VCF import VCFSmall
from VCFGeno import VCFGenoBase
from GenoRecord import GenoRecord
from VCFFamily import *
from VCFImpFamilyRecord import FillType, VCFImpFamilyRecord
from Genotype import Genotype
from TypeDeterminer import ParentComb
from Map import *
import Imputer
from option import *
from inverse_graph import *
from graph import Node
from invgraph import InvGraph


#################### VCFHeteroHomoRecord ####################

class VCFHeteroHomoRecord(VCFImpFamilyRecord):
	def __init__(self, pos: int, geno: list[int],
							i: int, parents_wrong_type: str,
							pair: ParentComb) -> None:
		super().__init__(pos, geno, i, parents_wrong_type, pair)
		self.which_comes_from: list[int] = [-1] * self.num_progenies()
	
	def is_mat_hetero(self) -> bool:
		return self.unphased_mat() == 1
	
	def is_pat_hetero(self) -> bool:
		return self.unphased_pat() == 1
	
	def is_imputable(self) -> bool:
		return self.parents_wrong_type == 'Right'
	
	def get_fill_type(self) -> FillType:
		if self.is_imputable():
			return FillType.MAT if self.is_mat_hetero() else FillType.PAT
		else:
			return FillType.IMPUTABLE
	
	def alleles_from_hetero_parent(self) -> list[int]:
		def encode(gt: int, homo_parent: int) -> int:
			gt_ = gt if homo_parent == Genotype.UN_00 else gt - 1
			if gt_ not in (0, 1):
				return -1
			else:
				return gt_
		
		gts = self.unphased_gts()
		# homoのGenotypeは0/0か1/1か
		homo = (self.unphased_pat() if self.is_mat_hetero()
									else self.unphased_mat())
		return [ encode(gt, homo) for gt in gts[2:] ]
	
	def set_haplo(self, h: int) -> None:
		hetero_index = 0 if self.is_mat_hetero() else 1
		homo_index = 1 if hetero_index == 0 else 0
		self.geno[hetero_index] = Genotype.PH_01 if h == 0 else Genotype.PH_10
		self.geno[homo_index] = (Genotype.PH_00
								 if self.unphased(homo_index) == Genotype.UN_00
								 else Genotype.PH_11)
		
		gts = self.alleles_from_hetero_parent()
		for i, gt in enumerate(gts):
			self.which_comes_from[i] = -1 if gt == -1 else 0 if gt == h else 1
	
	def set_int_gt_by_which_comes_from(self, ws: list[int]) -> None:
		# 親はphasingされている前提
		self.which_comes_from = ws
		mat_gt: int = self.mat_gt()
		pat_gt: int = self.pat_gt()
		if self.is_mat_hetero():
			pat_allele = pat_gt & 1
			for i in range(self.num_progenies()):
				mat_allele = (mat_gt >> ws[i]) & 1
				self.geno[i+2] = Genotype.from_alleles(mat_allele, pat_allele)
		else:
			mat_allele = mat_gt & 1
			for i in range(self.num_progenies()):
				pat_allele = (pat_gt >> ws[i]) & 1
				self.geno[i+2] = Genotype.from_alleles(mat_allele, pat_allele)


#################### VCFHeteroHomo ####################

class VCFHeteroHomo(VCFFamilyBase, VCFMeasurable):
	def __init__(self, samples: list[str], records: list[VCFHeteroHomoRecord],
											map_: Map, vcf: VCFSmall) -> None:
		VCFFamilyBase.__init__(self, samples, vcf)
		VCFMeasurable.__init__(self, map_)
		self.records: list[VCFHeteroHomoRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def is_mat_hetero(self) -> bool:
		return self.records[0].is_mat_hetero()
	
	def make_graph(self, max_dist: float) -> InvGraph:
		# ベータ分布を使うので、distはfloat
		def distance(gts1: list[int], gts2: list[int]) -> tuple[float, bool]:
			def dist_with_NA(right: int, counter_NA: int) -> float:
				diff = N - right - counter_NA
				diff_ratio = (diff + 1.0) / (right + diff + 2.0)
				return diff + counter_NA * diff_ratio
			
			N = len(gts1)
			counter_right = 0
			counter_diff = 0
			counter_NA = 0
			for int_gt1, int_gt2 in zip(gts1, gts2):
				if int_gt1 == -1 or int_gt2 == -1:
					counter_NA += 1
				elif int_gt1 == int_gt2:
					counter_right += 1
				else:
					counter_diff += 1
			
			if counter_right >= counter_diff:
				return (dist_with_NA(counter_right, counter_NA), False)
			else:
				return (dist_with_NA(counter_diff, counter_NA), True)
		
		L = len(self)
		
		# このグラフを作る処理は全てのレコード対で行うと時間がかかるので、
		# 10cM以内しかレコード対が繋がっているか調べない
		# しかし、レコードが少ないとその範囲でグラフが繋がるかわからないので、
		# 30レコードまでなら全てのレコード対を調べる
		# レコード数がそれ以上ならcMを考慮しなければ計算量が同じになるようにする
		# ただし、前後10レコード以上は調べる」
		rng = max(10, min(L, 900//L))
		
		graph: InvGraph = InvGraph()
		for k in range(L):
			graph[Node(k)] = []
		gtss = [ r.alleles_from_hetero_parent() for r in self.records ]
		cMs = [ self.cM(record.pos) for record in self.records ]
		# C++版ではkで並列化したい
		for k in range(L):
			for l in range(k+1, L):
				cM = cMs[l] - cMs[k]
				# 10cM以上離れていたら繋がりを見ない
				# ただし、一定の数のレコード分しか離れていなかったら見る
				if cM > 10.0 and k + rng < l:
					break
				d, b = distance(gtss[k], gtss[l])
				if d <= max_dist:
					graph[Node(k)].append((Node(l), d, b))
					graph[Node(l)].append((Node(k), d, b))
		return graph
	
	# haplotype1の各レコードのGenotypeが0なのか1なのか
	def make_parent_haplotypes(self, graph: InvGraph) -> list[int]:
		def is_visited(v: int) -> bool:
			return haplo[v] != -1
		
		tree: InvGraph = graph.minimum_spanning_tree()
		vs = sorted(tree.keys())
		haplo = self.create_haplotype(vs[0], tree)
		return [ h for h in haplo if h != -1 ]
	
	def create_haplotype(self, v0: Node, tree: InvGraph) -> list[int]:
		L = len(self)
		haplo: list[int] = [-1] * L
		haplo[v0] = 0
		edges = tree.walk(v0)
		for v1, v2, inv in edges:
			if haplo[v2] == -1:
				if inv:
					haplo[v2] = 1 if haplo[v1] == 0 else 0
				else:
					haplo[v2] = haplo[v1]
			elif haplo[v1] == -1:
				if inv:
					haplo[v1] = 1 if haplo[v2] == 0 else 0
				else:
					haplo[v1] = haplo[v2]
		return haplo
	
	def make_subvcf(self, graph: InvGraph) -> VCFHeteroHomo:
		haplo = self.make_parent_haplotypes(graph)
		indices = sorted(graph.keys())
		records = [ self.records[i] for i in indices ]
		for record, h in zip(records, haplo):
			record.set_haplo(h)
		
		return VCFHeteroHomo(self.samples, records, self.map, self.vcf)
	
	def determine_haplotype(self, option: OptionImpute
					) -> tuple[list[VCFHeteroHomo], list[VCFHeteroHomoRecord]]:
		max_dist = min(option.max_dist, (self.num_samples() - 2) * 0.1) * 2
		graph = self.make_graph(max_dist)
		subgraphs = graph.divide_into_connected()
		
		# 小さなgraphのmarkerは省く
		# 小さなgraphしかなければ、その中でも一番大きなgraphにする
		gs = [ g for g in subgraphs if len(g) >= option.min_graph ]
		if not gs:
			max_index, g = max(enumerate(subgraphs), key=lambda v: len(v[1]))
			gs = [g]
			# gに使われたrecordは入れない
			unused_records = [ self.records[j]
								for i, g in enumerate(subgraphs)
								if i != max_index and len(g) < option.min_graph
								for j in g.keys() ]
		else:
			unused_records = [ self.records[i]
							for g in subgraphs if len(g) < option.min_graph
							for i in g.keys() ]
		
		subvcfs = list(map(self.make_subvcf, gs))
		return (subvcfs, unused_records)
	
	def make_seq(self, i: int) -> str:
		def convert(gt: int) -> str:
			if gt == 0 or gt == 1:
				return str(gt)
			else:
				return 'N'
		
		return ''.join(convert(record.which_comes_from[i])
									for record in self.records)
	
	def impute_each_sample_seq(self, i: int, cMs: list[float],
													min_c: float) -> str:
		seq = self.make_seq(i)
		if Imputer.is_all_same_without_N(seq):
			return Imputer.create_same_color_string(seq, '0')
		
		hidden_states = ['0', '1']
		states = ['0', '1', 'N']
		hidden_seq = Imputer.impute(seq, hidden_states, states, cMs)
		painted_seq = Imputer.paint(hidden_seq, cMs, min_c)
		return painted_seq
	
	def impute_each(self, option: OptionImpute) -> None:
		cMs = [ self.cM(record.pos) for record in self.records ]
		imputed_seqs = [
				self.impute_each_sample_seq(i, cMs, option.min_crossover)
								for i in range(self.num_samples() - 2) ]
		for k, record in enumerate(self.records):
			ws = [ int(seq[k]) for seq in imputed_seqs ]
			record.set_int_gt_by_which_comes_from(ws)
	
	def create_option(self) -> OptionImpute:
		num = max(2, len(self))
		total_cM = self.cM(self.records[-1].pos)
		max_dist = max(4, int(total_cM * self.num_progenies()
											/ num * log10(num) * 2.5 * 0.01))
		return OptionImpute(max_dist, 20, 5, 1.0)
	
	def impute(self) -> tuple[list[VCFHeteroHomo], list[VCFHeteroHomoRecord]]:
		if not self.records:
			return ([self], [])
		option: OptionImpute = self.create_option()
		vcfs, unused_records = self.determine_haplotype(option)
		for vcf in vcfs:
			vcf.impute_each(option)
		for record in unused_records:
			record.enable_modification()
		return (vcfs, unused_records)
	
	# 共通のヘテロ親はどれだけマッチしているか
	def match(self, other: VCFHeteroHomo) -> tuple[int, int]:
		if len(self) == 0 or len(other) == 0:
			return (0, 0)
		
		hetero_col1 = 0 if self.is_mat_hetero() else 1
		hetero_col2 = 0 if other.is_mat_hetero() else 1
		L, M = len(self), len(other)
		num_match, num_unmatch = 0, 0
		k, l = 0, 0
		while k < L and l < M:
			record1 = self.records[k]
			record2 = other.records[l]
			if record1.pos == record2.pos:
				if record1.geno[hetero_col1] == record2.geno[hetero_col2]:
					num_match += 1
				else:
					num_unmatch += 1
			if record1.pos <= record2.pos:
				k += 1
			if record1.pos >= record2.pos:
				l += 1
		
		return (num_match, num_unmatch)
	
	def inverse_hetero_parent_phases(self) -> None:
		hetero_col = 0 if self.is_mat_hetero() else 1
		for record in self.records:
			if record.geno[hetero_col] == Genotype.PH_01:
				record.geno[hetero_col] = Genotype.PH_10
			else:
				record.geno[hetero_col] = Genotype.PH_01
	
	# Create a pair of VCFHeteroHomo for each Family and store it for each parent
	@staticmethod
	def make_VCFHeteroHomo(records: list[VCFHeteroHomoRecord],
							samples: list[str], geno_map: Map, vcf: VCFSmall
									) -> tuple[VCFHeteroHomo, VCFHeteroHomo,
													list[VCFHeteroHomoRecord]]:
		heho_mat_records = [ r for r in records
								if r.is_imputable() and not r.is_pat_hetero() ]
		heho_pat_records = [ r for r in records
								if r.is_imputable() and not r.is_mat_hetero() ]
		# TODO: classifyのときにunusedに入れておけばいいのでは？
		unused_records = [ r for r in records if not r.is_imputable() ]
		vcf_mat = VCFHeteroHomo(samples, heho_mat_records, geno_map, vcf)
		vcf_pat = VCFHeteroHomo(samples, heho_pat_records, geno_map, vcf)
		return (vcf_mat, vcf_pat, unused_records)
	
	@staticmethod
	def impute_vcfs(records: list[VCFHeteroHomoRecord],
					samples: list[str], geno_map: Map, vcf: VCFSmall
					) -> tuple[list[VCFHeteroHomo], list[VCFHeteroHomoRecord]]:
		vcf_mat, vcf_pat, unused = VCFHeteroHomo.make_VCFHeteroHomo(
											records, samples, geno_map, vcf)
		vcfs_mat, unused_records1 = vcf_mat.impute()
		vcfs_pat, unused_records2 = vcf_pat.impute()
		vcfs = []
		vcfs.extend(vcfs_mat)
		vcfs.extend(vcfs_pat)
		unused.extend(unused_records1)
		unused.extend(unused_records2)
		return (vcfs, unused)
	
	@staticmethod
	def inverse_phases(vcfs: list[VCFHeteroHomo]) -> None:
		# 同じ親が複数の家系にまたがるとき
		graph = VCFHeteroHomo.make_vcf_graph(vcfs)
		if not graph:
			return
		
		bs: tuple[bool, ...] = graph.optimize_inversions()
		for vcf, b in zip(vcfs, bs):
			if b:
				vcf.inverse_hetero_parent_phases()
	
	@staticmethod
	def make_vcf_graph(vcfs: list[VCFHeteroHomo]) -> InverseGraph:
		N = len(vcfs)
		graph: InverseGraph = InverseGraph()
		for v in range(N):
			graph[Node(v)] = []
		for i, j in combinations(range(N), 2):
			num_match, num_unmatch = vcfs[i].match(vcfs[j])
			if num_match != 0 or num_unmatch != 0:
				graph[Node(i)].append((Node(j), num_match, num_unmatch))
				graph[Node(j)].append((Node(i), num_match, num_unmatch))
		return graph

__all__ = ['VCFHeteroHomoRecord', 'VCFHeteroHomo']
