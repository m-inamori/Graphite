from __future__ import annotations

# coding: utf-8
# VCFSelfHetero.py
# 自殖で親がヘテロ

from collections import defaultdict, Counter
from math import log10
from typing import List, Tuple, Optional, IO, Dict, Iterator
import random
import time

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from GenoRecord import GenoRecord
from VCFImpSelfRecord import SelfFillType, VCFImpSelfRecord
from VCFSelfHeteroRecord import VCFSelfHeteroRecord
from SelfProgenyImputer import SelfProgenyImputer
from TypeDeterminer import ParentComb
from Map import *
import Imputer
from Genotype import Genotype
from option import *
from inverse_graph import *
from graph import Node
from invgraph import InvGraph


#################### VCFSelfHetero ####################

class VCFSelfHetero(VCFGenoBase, VCFMeasurable):
	def __init__(self, samples: list[str], records: list[VCFSelfHeteroRecord],
											map_: Map, vcf: VCFSmall) -> None:
		VCFGenoBase.__init__(self, samples, vcf)
		VCFMeasurable.__init__(self, map_)
		self.records = records
		self.prog_imputer = SelfProgenyImputer(records, map_, 0.01)
	
	def __len__(self) -> int:
		return len(self.records)
	
	def num_progenies(self) -> int:
		return self.num_samples() - 1
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	def make_graph(self, max_dist: float) -> InvGraph:
		# ベータ分布を使うので、distはfloat
		def distance(gts1: list[int], gts2: list[int]) -> tuple[float, bool]:
			def dist_with_NA(right: int, counter_NA: int) -> float:
				diff = N - right - counter_NA
				diff_ratio = (diff + 1.0) / (right + diff + 2.0)
				return diff + counter_NA * diff_ratio
			
			N = len(gts1)
			counter_homo_right = 0
			counter_hetero_right = 0
			counter_inverse = 0
			counter_diff = 0
			counter_NA = 0
			for gt1, gt2 in zip(gts1, gts2):
				if gt1 == Genotype.NA or gt2 == Genotype.NA:
					counter_NA += 1
				elif (gt1, gt2) == (Genotype.UN_01, Genotype.UN_01):
					counter_hetero_right += 1
				elif gt1 == gt2:
					counter_homo_right += 1
				elif (gt1, gt2) in ((0, 2), (2, 0)):
					counter_inverse += 1
				else:
					counter_diff += 1
			
			if counter_homo_right >= counter_inverse:
				right = counter_hetero_right + counter_homo_right
				return (dist_with_NA(right, counter_NA), False)
			else:
				right = counter_hetero_right + counter_inverse
				return (dist_with_NA(right, counter_NA), True)
		
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
		gtss = [ [ Genotype.all_int_gt_to_int_gt(r.geno[i])
							 for i in range(1, self.num_samples()) ]
													 for r in self.records ]
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
	
	def make_subvcf(self, graph: InvGraph) -> VCFSelfHetero:
		haplo = self.make_parent_haplotypes(graph)
		indices = sorted(graph.keys())
		records = [ self.records[i] for i in indices ]
		for record, h in zip(records, haplo):
			record.set_haplo(h)
		
		return VCFSelfHetero(self.samples, records, self.map, self.vcf)
	
	def determine_haplotype(self, option: OptionImpute
										) -> tuple[list[VCFSelfHetero],
												   list[VCFSelfHeteroRecord]]:
		max_dist = min(option.max_dist, (self.num_samples() - 1) * 0.1) * 2
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
	
	def impute_progenies(self, option: OptionImpute) -> None:
		for iprog in range(self.num_progenies()):
			self.prog_imputer.impute(iprog)
	
	def create_option(self) -> OptionImpute:
		num = max(2, len(self))
		total_cM = self.cM(self.records[-1].pos)
		max_dist = max(4, int(total_cM * self.num_progenies()
											/ num * log10(num) * 2.5 * 0.01))
		return OptionImpute(max_dist, 20, 5, 1.0)
	
	def impute(self) -> tuple[list[VCFSelfHetero], list[VCFSelfHeteroRecord]]:
		if not self.records:
			return ([self], [])
		option: OptionImpute = self.create_option()
		vcfs, unused_records = self.determine_haplotype(option)
		for vcf in vcfs:
			vcf.impute_progenies(option)
		for record in unused_records:
			record.enable_modification()
		return (vcfs, unused_records)

__all__ = ['VCFSelfHeteroRecord', 'VCFSelfHetero']
