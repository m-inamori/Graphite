from __future__ import annotations

# coding: utf-8
# VCFOneParentImputedMCMC.py
# 片側だけimputeされている
# MCMCを使う

from collections import defaultdict, Counter
from math import log
from typing import Tuple, Dict, Iterator
import copy
import random

from VCFFamily import *
from Map import *
from Genotype import Genotype


#################### State ####################

class State:
	def __init__(self, s: list[list[int]], p: float):
		self.haps = s
		self.prob = p


#################### VCFOneParentImputedMCMC ####################

class VCFOneParentImputedMCMC(VCFBase, VCFSmallBase,
								VCFFamilyBase, VCFMeasurable):
	def __init__(self, header: list[list[str]], records: list[VCFFamilyRecord],
					ref_haps: list[list[int]], is_mat_imputed: bool, map_: Map):
		VCFBase.__init__(self, header)
		VCFSmallBase.__init__(self)
		VCFFamilyBase.__init__(self)
		VCFMeasurable.__init__(self, map_)
		self.records: list[VCFFamilyRecord] = records
		self.ref_haps = ref_haps
		self.is_mat_imputed = is_mat_imputed
		self.phased_col = 9 if self.is_mat_imputed else 10
		self.non_phased_col = 10 if self.is_mat_imputed else 9
		# crossover values
		# 後代
		self.Cc = [ Map.Kosambi(self.dist(r1, r2))
							for r1, r2 in zip(self.records, self.records[1:]) ]
		K = 5		# 親の遷移確率の子に対する倍率
		# 親は直接の後代でないので、乗り換え確率を高くする
		self.Cp = [ Map.Kosambi(self.dist(r1, r2) * K)
							for r1, r2 in zip(self.records, self.records[1:]) ]
		self.phased_parent_gts = [ Genotype.phased_gt_to_int(
												r.v[self.phased_col])
														for r in self.records ]
		self.observed_parent = [ Genotype.gt_to_int(r.v[self.non_phased_col])
														for r in self.records ]
		self.observed_progs = [ [ Genotype.gt_to_int(gt) for gt in r.v[11:] ]
														for r in self.records ]
		self.compute_exhaust_probabilities()
	
	def compute_exhaust_probabilities(self):
		# 排出確率
		# hidden   0|0: 0, 1|0: 1, 0|1: 2, 1|1: 3
		# observed 0/0: 0, 0/1: 1, 1/1: 2, N/A: 3
		w = 0.01		# probability of wrong
		E1 = [ [ 1.0-3*w, w,       w,       w ],
			   [ w,       1.0-3*w, w,       w ],
			   [ w,       1.0-3*w, w,       w ],
			   [ w,       w,       1.0-3*w, w ] ]
		
		self.E = [ [ log(p) for p in v ] for v in E1 ]
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> VCFRecord:
		return self.records[i]
	
	def get_family_record(self, i: int) -> VCFFamilyRecord:
		return self.records[i]
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def dist(self, r1: VCFRecord, r2: VCFRecord) -> float:
		dist = (self.cM(r2.pos()) - self.cM(r1.pos())) / 100
		if dist > 0.0:
			return dist
		else:
			return (r2.pos() - r1.pos()) * 1e-6
	
	def initialize_state(self) -> State:
		N = self.num_progenies()
		M = len(self)				# マーカー数
		haps = [ [0] * M for _ in range(N * 2 + 2) ]
		p = self.compute_probability(haps)
		return State(haps, p)
	
	def get_non_phased_parent_gt(self, hs: list[list[int]], i: int) -> int:
		hp1 = hs[0][i]
		hp2 = hs[1][i]
		return self.ref_haps[hp1][i] | (self.ref_haps[hp2][i] << 1)
	
	def emission_probability(self, hs: list[list[int]], i: int) -> float:
		record = self.records[i]
		phased_parent_gt = Genotype.phased_gt_to_int(record.v[self.phased_col])
		non_phased_parent_gt = self.get_non_phased_parent_gt(hs, i)
		mat_gt = (phased_parent_gt if self.is_mat_imputed
										else non_phased_parent_gt)
		pat_gt = (non_phased_parent_gt if self.is_mat_imputed
										else phased_parent_gt)
		Ep = self.E[non_phased_parent_gt][self.observed_parent[i]]
		# 後代の排出確率を計算
		Ec = 0.0
		for k in range(2, len(hs), 2):
			j = (k - 2) >> 1
			hc1 = hs[k][i]
			hc2 = hs[k+1][i]
			gtc = Genotype.gt_by_haplotypes(hc1, hc2, mat_gt, pat_gt)
			Ec += self.E[gtc][self.observed_progs[i][j]]
		return Ep + Ec
	
	def parent_transition_probability(self, hs: list[list[int]]) -> float:
		prob = 0.0
		for i in range(2):
			for j in range(len(hs[0]) - 1):
				if hs[i][j] != hs[i][j+1]:
					prob += log(self.Cp[j])
				else:
					prob += log(1.0 - self.Cp[j])
		return prob
	
	def progeny_transition_probability(self, hs: list[list[int]]) -> float:
		prob = 0.0
		for i in range(2, len(hs)):
			for j in range(len(hs[0]) - 1):
				if hs[i][j] != hs[i][j+1]:
					# 前もってlogの計算をしたい
					prob += log(self.Cc[j])
				else:
					prob += log(1.0 - self.Cc[j])
		return prob
	
	def transition_probability(self, hs: list[list[int]]) -> float:
		return (self.parent_transition_probability(hs) +
				self.progeny_transition_probability(hs))
	
	def compute_probability(self, hs: list[list[int]]) -> float:
		# 計算量低減のため、変えた部分のみ再計算するようにしたい
		return (self.transition_probability(hs) + 
				sum(self.emission_probability(hs, i) for i in range(len(self))))
	
	def join_randomly(self, state: State) -> State:
		M = len(self)
		divided_points = [ (i, j) for i, hap in enumerate(state.haps)
								  for j in range(M-1) if hap[j] != hap[j+1] ]
		if not divided_points:
			return state	# ジョインできないのでそのまま
		
		i, j = random.choice(divided_points)
		b = random.randrange(0, 2) == 0
		haps = copy.copy(state.haps)
		if b:
			# 前半と同じにする
			c = haps[i][j]
			c1 = haps[i][j+1]
			for j1 in range(j+1, M):
				if haps[i][j1] != c1:
					break
				haps[i][j1] = c
		else:
			c = state.haps[i][j+1]
			c1 = state.haps[i][j]
			for j1 in range(j, -1, -1):
				if haps[i][j1] != c1:
					break
				haps[i][j1] = c
		p = self.compute_probability(haps)
		return State(haps, p)
	
	def nearest_crossover_point_lower(self, s: State, i: int, j: int) -> int:
		for j1 in range(j + 1, len(self) - 1):
			if s.haps[i][j1] != s.haps[i][j1+1]:
				return j1
		return len(self)
	
	def select_different_color(self, haps: list[list[int]],
												i: int, j: int) -> int:
		c0 = haps[i][j]
		if i > 1:
			return 1 if c0 == 0 else 0
		else:
			c = random.randrange(len(self.ref_haps) - 1)
			return c if c < c0 else c + 1
	
	def change_upper_color(self, c: int, i: int, j: int, haps: list[list[int]]):
		c0 = haps[i][j]
		for j1 in range(j - 1, -1, -1):
			if haps[i][j1] != c0:
				break
			haps[i][j1] = c
	
	def change_lower_color(self, c: int, i: int, j: int, haps: list[list[int]]):
		M = len(self)
		c0 = haps[i][j]
		for j1 in range(j + 1, M):
			if haps[i][j1] != c0:
				break
			haps[i][j1] = c
	
	def decide_new_division_point(self, state: State) -> tuple[int, int]:
		M = len(self)
		while True:
			N = len(state.haps)
			r = random.randrange((M - 1) * N)
			i, j = divmod(r, M - 1)		# i番目のhapのj番目のマーカーで分割
			if state.haps[i][j] == state.haps[i][j+1]:
				return (i, j)
	
	def divide_randomly(self, state: State) -> State:
		# 浅いコピーにして、変えるハプロタイプのみ深いコピーにする
		haps = copy.copy(state.haps)
		i, j = self.decide_new_division_point(state)
		haps[i] = copy.copy(haps[i])
		c = self.select_different_color(haps, i, j)
		if random.randrange(0, 2) == 0:
			self.change_upper_color(c, i, j, haps)	# 上の色を変える
			haps[i][j] = c
		else:
			self.change_lower_color(c, i, j, haps)	# 下の色を変える
		p = self.compute_probability(haps)
		return State(haps, p)
	
	def change_color_randomly(self, state: State) -> State:
		haps = copy.copy(state.haps)
		M = len(self)
		N = len(state.haps)
		r = random.randrange(M * N)
		# i番目のhapのj番目のマーカーを含む領域の色を変える
		i, j = divmod(r, M)
		haps[i] = copy.copy(haps[i])
		c = self.select_different_color(haps, i, j)
		self.change_upper_color(c, i, j, haps)	# 上の色を変える
		self.change_lower_color(c, i, j, haps)	# 下の色を変える
		haps[i][j] = c
		p = self.compute_probability(haps)
		return State(haps, p)
	
	def proposal_func(self, state: State) -> State:
		r = random.random()
		if r < 1 / 3:
			return self.join_randomly(state)
		elif r < 2 / 3:
			return self.divide_randomly(state)
		else:
			return self.change_color_randomly(state)
	
	def mcmc_sampling_hmm(self, initial_state: State, num_iters: int) -> State:
		current_state = initial_state
		best_state = current_state
		
		for _ in range(num_iters):
			proposed_state = self.proposal_func(current_state)
			ratio = proposed_state.prob - current_state.prob
			if log(random.random()) < ratio:
				current_state = proposed_state
				if proposed_state.prob < best_state.prob:
					best_state = current_state
		
		return best_state
	
	def compute_non_phased_parent_gt(self, hs: list[list[int]], i: int) -> int:
		hp1 = hs[0][i]
		hp2 = hs[1][i]
		return self.ref_haps[hp1][i] | (self.ref_haps[hp2][i] << 1)
	
	def update_genotypes(self, state: State):
		for i, record in enumerate(self.records):
			hs = state.haps
			non_phased_parent_gt = self.compute_non_phased_parent_gt(hs, i)
			record.v[self.non_phased_col] = Genotype.int_to_phased_gt(
														non_phased_parent_gt)
			gt = record.v[self.phased_col]
			phased_parent_gt = Genotype.phased_gt_to_int(gt)
			mat_gt = (phased_parent_gt if self.is_mat_imputed
											else non_phased_parent_gt)
			pat_gt = (non_phased_parent_gt if self.is_mat_imputed
											else phased_parent_gt)
			for j in range(2, len(hs), 2):	# 個々の後代
				hc1 = hs[j][i]
				hc2 = hs[j+1][i]
				gtc_int = Genotype.gt_by_haplotypes(hc1, hc2, mat_gt, pat_gt)
				record.v[j//2+10] = Genotype.int_to_phased_gt(gtc_int)
	
	def impute(self):
		random.seed(2)
		initial_state = self.initialize_state()
		state = self.mcmc_sampling_hmm(initial_state, 100000)
		self.update_genotypes(state)
