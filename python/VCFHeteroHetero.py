from __future__ import annotations
from collections import defaultdict

from VCFFamily import *
from Map import *


#################### VCFHeteroHeteroRecord ####################

class VCFHeteroHeteroRecord(VCFFamilyRecord):
	def __init__(self, v: list[str], samples: list[str], index, w1=[], w2=[]):
		super().__init__(v, samples)
		self.index: int = index
		self.which_comes_from1: list[str] = w1
		self.which_comes_from2: list[str] = w2
	
	def update_genotype(self, i: int):
		if self.which_comes_from1[i] == '2':
			gt = '0/1'
		elif self.v[9] == '0|1':
			gt = self.which_comes_from1[i] + '|' + self.which_comes_from2[i]
		else:
			g1 = '0' if self.which_comes_from1[i] == '1' else '1'
			g2 = '0' if self.which_comes_from2[i] == '1' else '1'
			gt = g1 + '|' + g2
		self.v[i+11] = gt + self.v[i+11][3:]
	
	# (0|1, 0|1), mat_gt, pat_gt -> gt
	@staticmethod
	def make_genotype(which: tuple[int,int], F1_gt: str) -> str:
		w1, w2 = which
		return F1_gt[w1*2] + '|' + F1_gt[w2*2]


#################### VCFHeteroHetero ####################

class VCFHeteroHetero(VCFFamilyBase, VCFMeasurable):
	class OptionImpute:
		def __init__(self, max_dist, min_positions, min_graph):
			self.max_dist: int		= max_dist
			self.min_positions: int	= min_positions
			self.min_graph: int		= min_graph
	
	def __init__(self, header: list[list[str]],
					records: list[VCFHeteroHeteroRecord], map_: Map):
		self.records = records
		VCFFamilyBase.__init__(self, header)
		VCFMeasurable.__init__(self, map_)
	
	# [int], [int], float -> (dist: float, change type: (0|1|2|-1))
	# ベータ分布を使うので、distはfloat
	def distance(self, gts1: list[int], gts2: list[int],
									max_dist: float) -> tuple[float, int]:
		def dist_with_NA(right: int, counter_NA: int) -> float:
			diff = N - right - counter_NA
			diff_ratio = (diff + 1.0) / (right + diff + 2.0)
			return diff + counter_NA * diff_ratio
		
		N = len(gts1)
		# int_gt1 -> int_gt2の頻度
		counter: list[int] = [0] * 9	# 0->0, 0->1, ... 2->2
		counter_NA = 0
		for int_gt1, int_gt2 in zip(gts1, gts2):
			if int_gt1 == -1 or int_gt2 == -1:
				counter_NA += 1
			else:
				counter[int_gt1+int_gt2*3] += 1
		
		# type 0 (not change)
		right = counter[0] + counter[4] + counter[8]
		d0: float = dist_with_NA(right, counter_NA)
		if d0 <= max_dist:
			return (d0, 0)
		
		# type 1 (0 -> 1, 0 <- 1 -> 2, 2 -> 1)
		right = counter[1] + counter[3] + counter[5] + counter[7]
		d1: float = dist_with_NA(right, counter_NA)
		if d1 <= max_dist:
			return (d1, 1)
		
		# type 2 (0 -> 2, 1 -> 1, 2 -> 0)
		right = counter[2] + counter[4] + counter[6]
		d2: float = dist_with_NA(right, counter_NA)
		if d2 <= max_dist:
			return (d2, 2)
		
		return (max_dist + 1.0, -1)
	
	def make_graph(self, option: OptionImpute
							) -> Dict[int, tuple[int, float, int]]:
		graph = defaultdict(list)
		L = len(self.records)
		cMs = [ self.map.bp_to_cM(self.records[i].pos()) for i in range(L) ]
		max_dist: float = min(option.max_dist, self.num_progenies() * 0.1)
		gtss = [ record.progeny_int_gts() for record in self.records ]
		for k in range(L):
			for l in range(k+1, L):
				cM = cMs[l] - cMs[k]
				if cM > 10.0:	# 10cM以上離れていたら繋がりを見ない
					break
				d, t = self.distance(gtss[k], gtss[l], max_dist)
				if d <= max_dist:
					graph[k].append((l, d, t))
					graph[l].append((k, d, t))
		return graph
	
	def determine_haplotype(self, graph):
		# minimum spanning tree
		def make_tree(graph):
			# { node: [(node, weight, type)] } -> { node: [(node, weight)] }
			weighted_graph = { v: [ t[:2] for t in vs ]
										for v, vs in graph.items() }
			tree = Kluskal_fast(weighted_graph)		# { node: [(node, weight)] }
			new_tree = { }
			for v, vs in tree.items():
				types = dict((v1, t) for v1, d, t in graph[v])
				new_tree[v] = [ (v1, types[v1]) for v1, d in vs ]
			return new_tree
		
		# v: int -> { node: prev node }
		def walk(tree, v0):
			visited = set([v0])
			stk = [v0]
			prevs = { -1: -1, v0: -1 }
			while stk:
				v = stk.pop()
				for elem in tree[v]:
					v1 = elem[0]
					if v1 in visited:
						continue
					stk.append(v1)
					prevs[v1] = v
					visited.add(v1)
			
			return prevs
		
		def modify_record(record, t):
			v = record.v[:]
			gt = '0|1' if t == 0 else '1|0'
			v[9] = gt + v[9][3:]
			v[10] = gt + v[10][3:]
			which1 = []
			which2 = []
			for c in range(11, len(v)):
				if record.v[c] == '0/0':
					v[c] = '0|0' + v[c][3:]
					w = 0 if gt[0] == '0' else 1
					which1.append(w)
					which2.append(w)
				elif record.v[c] == '1/1':
					v[c] = '1|1' + v[c][3:]
					w = 1 if gt[0] == '0' else 0
					which1.append(w)
					which2.append(w)
				elif record.v[c] == '0/1':
					which1.append(2)
					which2.append(2)
				else:				# ./.
					which1.append(-1)
					which2.append(-1)
			return VCFHeteroHeteroRecord(v, record.samples,
											record.index, which1, which2)
		
		def determine_parents(tree, v0):
			records = [None] * len(self)
			stk = [(v0, 0)]
			visited = set([v0])
			while stk:
				v, t = stk.pop()
				records[v] = modify_record(self.records[v], t)
				for v1, type in tree[v]:
					if v1 in visited:
						continue
					stk.append((v1, t ^ type))
					visited.add(v1)
			return records
		
		tree = make_tree(graph)
		v0 = min(tree.keys())
		return determine_parents(tree, v0)
	
	def make_subvcf(self, graph):
		records = self.determine_haplotype(graph)
		return VCFHeteroHetero(self.header, records, self.map)
	
	def impute_each_seq(self, seq, hidden_states, states, Ts):
		pi, E = BWTL.Baum_Welch(Ts, seq, states, hidden_states)
		A = BWTL.reverse_probs(E, hidden_states, states)
		hidden_seq = BWTL.Viterbi(seq, states, hidden_states, pi, Ts, A)
		return hidden_seq
#		painted_seq = self.paint(hidden_seq)
#		return painted_seq
	
	def make_seq(self, which):
		def convert(gt):
			if gt == -1:
				return 'N'
			else:
				return str(gt)
		
		return ''.join(map(convert, which))
	
	def impute_each(self):
		def compute_T(p):
			# 0 <-> 1は両方の染色体で乗り換えが起こるから確率は小さい
			q = 1.0 - p
			return { ('0', '0'): q**2, ('0', '1'): p*p, ('0', '2'): 2*p*q,
					 ('1', '0'): p*p, ('1', '1'): q**2, ('1', '2'): 2*p*q,
					 ('2', '0'): p*q, ('2', '1'): p*q, ('2', '2'): q**2 + p**2 }
		
		records = [ record for record in self.records if record is not None ]
		hidden_states = ['0', '1', '2']
		v = [ self.cM(i) for i in range(len(self))
								if self.records[i] is not None ]
		ps = [ (p2-p1)/100 for p1, p2 in zip(v, v[1:]) ]
		Ts = [ compute_T(p) for p in ps ]
		states = [ s for s in ['0', '1', '2', 'N'] ]
		for i in range(self.num_progenies()):
			which = [ record.which_comes_from1[i] for record in records ]
			seq = self.make_seq(which)
			imputed = self.impute_each_seq(seq, hidden_states, states, Ts)
			for record, w in zip(records, imputed):
				record.which_comes_from1[i] = w
				record.which_comes_from2[i] = w
				record.update_genotype(i)
	
	def impute(self):
		option = VCFHeteroHetero.OptionImpute(4, 20, 5)
		
		##### debug code #####
		
#		self.records = self.records[:1000]
		graph = self.make_graph(option)
		subgraphs = divide_graph_into_connected(graph)
		imputed_vcfs = []
		for subg in subgraphs:
			if len(subg) >= option.min_positions:
				subvcf = self.make_subvcf(subg)
				subvcf.impute_each()
				imputed_vcfs.append(subvcf)
		return imputed_vcfs
