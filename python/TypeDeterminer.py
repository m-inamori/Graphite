from __future__ import annotations

# coding: utf-8
# TypeDeterminer.py

from math import log
from itertools import product, count
from collections import defaultdict
from queue import PriorityQueue
from enum import Enum
from typing import Dict, Iterator

from MathCommon import Simpson


#################### ParentComb ####################

class ParentComb(Enum):
	P00x00 = 0
	P00x01 = 1
	P01x01 = 2
	P00x11 = 3
	P01x11 = 4
	P11x11 = 5
	PNA = 6
	
	def is_NA(self) -> bool:
		return self == ParentComb.PNA
	
	def is_homohomo(self) -> bool:
		return self in (ParentComb.P00x00, ParentComb.P00x11, ParentComb.P11x11)
	
	def is_heterohomo(self) -> bool:
		return self in (ParentComb.P00x01, ParentComb.P01x11)
	
	def is_same_parent_genotype(self) -> bool:
		return self in (ParentComb.P00x00, ParentComb.P01x01, ParentComb.P11x11)
	
	def int_gt_pair(self) -> tuple[int, int]:
		p = self.value
		if p == 0:
			return (0, 0)
		elif p < 3:
			return (p - 1, 1)
		else:
			return (p - 3, 2)


#################### TypeDeterminer ####################

class TypeDeterminer:
	def __init__(self, n: int):
		self.N: int = n
		self.a = 1	# エラーの事前分布のparameter
		self.b = 9
		self.log_facs = self.make_memo_log_facs()
	
	def make_memo_log_facs(self) -> list[float]:
		L = self.N + self.b + 3
		v = [0.0] * L
		for n in range(1, L):
			v[n] = v[n-1] + log(n)
		return v
	
	def log_beta(self, n: int, m: int) -> float:
		return self.log_facs[n-1] + self.log_facs[m-1] - self.log_facs[n+m-1]
	
	# エラーが無ければ全て同じ遺伝子型になる場合
	def log_likelihood_one(self, n: int, m: int) -> float:
		return (self.log_beta(m+self.a, n+self.b) -
				m*log(3.0) - self.log_beta(self.a, self.b))
	
	# log likelihood for 0/0 x 0/0
	def log_likelihood00(self, mat_gt: int, pat_gt: int,
										gt_freq: list[int]) -> float:
		# 0/0 x 0/0 -> 0/0: 1-p otherwise: p/3
		n = gt_freq[0]								# 0/0
		m = gt_freq[1] + gt_freq[2] + gt_freq[3]	# 0/1 + 1/1 + ./.
		pair_gts = (mat_gt, pat_gt)
		if pair_gts == (0, 0):
			# 両親合っている
			# 両親の尤度は(1-p)^2
			return self.log_likelihood_one(n+2, m)
		elif mat_gt == 0 or pat_gt == 0:
			# 片親だけ合っている
			# 両親の尤度は(1-p)p
			return self.log_likelihood_one(n+1, m+1)
		else:
			# 両親とも合っていない
			# 両親の尤度はp^2
			return self.log_likelihood_one(n, m+2)
	
	# log likelihood for 0/0 x 0/1
	def log_likelihood01(self, mat_gt: int, pat_gt: int,
										gt_freq: list[int]) -> float:
		# 0/0 x 0/1 ->	-> 0/0 or 0/1: 1/2-p/3, otherwise: p/3
		n = gt_freq[0] + gt_freq[1]		# 0/0 or 0/1
		m = gt_freq[2] + gt_freq[3]		# 1/1 or ./.
		pair_gts = (mat_gt, pat_gt)
		if pair_gts == (0, 1) or pair_gts == (1, 0):
			# 両親合っている
			# 両親の尤度は(1-p)^2
			f = lambda p: (1-p)**2 * (1/2-p/3)**n * (p/3)**m * (1-p)**8
			return log(Simpson(f, 0.0, 1.0, 10)) - self.log_beta(self.a, self.b)
		elif mat_gt == 0 or mat_gt == 1 or pat_gt == 0 or pat_gt == 1:
			# 片親だけ合っている
			# 両親の尤度は(1-p)p/3
			f = lambda p: (1-p)*p/3 * (1/2-p/3)**n * (p/3)**m * (1-p)**8
			return log(Simpson(f, 0.0, 1.0, 10)) - self.log_beta(self.a, self.b)
		else:
			# 両親とも合っていない
			# 両親の尤度は(p/3)^2
			f = lambda p: (p/3)**2 * (1/2-p/3)**n * (p/3)**m * (1-p)**8
			return log(Simpson(f, 0.0, 1.0, 10)) - self.log_beta(self.a, self.b)
	
	# log likelihood for 0/1 x 0/1
	def log_likelihood11(self, mat_gt: int, pat_gt: int,
										gt_freq: list[int]) -> float:
		# 0/1 x 0/1 -> 0/0: 1/4, 0/1: 1/2-p/3, 1/1: 1/4, otherwise: p/3
		n = gt_freq[0] + gt_freq[2]
		m = gt_freq[1]
		q = gt_freq[3]
		pair_gts = (mat_gt, pat_gt)
		if pair_gts == (1, 1):
			# 両親合っている
			# 両親の尤度は(1-p)^2
			f = lambda p: (1-p)**2 * (1/4)**n * (1/2-p/3)**m * (p/3)**q * (1-p)**8
			return log(Simpson(f, 0.0, 1.0, 10)) - self.log_beta(self.a, self.b)
		elif mat_gt == 1 or pat_gt == 1:
			# 片親だけ合っている
			# 両親の尤度は(1-p)p/3
			f = lambda p: (1-p)*p/3 * (1/4)**n * (1/2-p/3)**m * (p/3)**q * (1-p)**8
			return log(Simpson(f, 0.0, 1.0, 10)) - self.log_beta(self.a, self.b)
		else:
			# 両親とも合っていない
			# 両親の尤度は(p/3)^2
			f = lambda p: (p/3)**2 * (1/4)**n * (1/2-p/3)**m * (p/3)**q * (1-p)**8
			return log(Simpson(f, 0.0, 1.0, 10)) - self.log_beta(self.a, self.b)
	
	# log likelihood for 0/0 x 1/1
	def log_likelihood02(self, mat_gt: int, pat_gt: int,
										gt_freq: list[int]) -> float:
		# 0/0 x 1/1 ->	-> 0/1: 1-p, otherwise: p/3
		n = gt_freq[1]
		m = gt_freq[0] + gt_freq[2] + gt_freq[3]
		pair_gts = (mat_gt, pat_gt)
		if pair_gts == (0, 2) or pair_gts == (2, 0):
			# 両親合っている
			# 両親の尤度は(1-p)^2
			return self.log_likelihood_one(n+2, m)
		elif mat_gt == 0 or mat_gt == 2 or pat_gt == 0 or pat_gt == 2:
			# 片親だけ合っている
			# 両親の尤度は(1-p)p/3
			return self.log_likelihood_one(n+1, m+1)
		else:
			# 両親とも合っていない
			# 両親の尤度は(p/3)^2
			return self.log_likelihood_one(n, m+2)
	
	# log likelihood for 0/1 x 1/1
	def log_likelihood12(self, mat_gt: int, pat_gt: int,
										gt_freq: list[int]) -> float:
		# 0/1 x 1/1 -> 0/1 or 1/1: 1/2-p/3, otherwise: p/3
		n = gt_freq[1] + gt_freq[2]		# 0/1 or 1/1
		m = gt_freq[0] + gt_freq[3]		# 0/0 or ./.
		pair_gts = (mat_gt, pat_gt)
		if pair_gts == (1, 2) or pair_gts == (2, 1):
			# 両親合っている
			# 両親の尤度は(1-p)^2
			f = lambda p: (1-p)**2 * (1/2-p/3)**n * (p/3)**m * (1-p)**8
			return log(Simpson(f, 0.0, 1.0, 10)) - self.log_beta(self.a, self.b)
		elif mat_gt == 1 or mat_gt == 2 or pat_gt == 1 or pat_gt == 2:
			# 片親だけ合っている
			# 両親の尤度は(1-p)p/3
			f = lambda p: (1-p)*p/3 * (1/2-p/3)**n * (p/3)**m * (1-p)**8
			return log(Simpson(f, 0.0, 1.0, 10)) - self.log_beta(self.a, self.b)
		else:
			# 両親とも合っていない
			# 両親の尤度はp^2
			f = lambda p: (p/3)**2 * (1/2-p/3)**n * (p/3)**m * (1-p)**8
			return log(Simpson(f, 0.0, 1.0, 10)) - self.log_beta(self.a, self.b)
	
	# log likelihood for 1/1 x 1/1
	def log_likelihood22(self, mat_gt: int, pat_gt: int,
										gt_freq: list[int]) -> float:
		# 1/1 x 1/1 -> 1/1: 1-p otherwise: p/3
		n = gt_freq[2]								# 1/1
		m = gt_freq[0] + gt_freq[1] + gt_freq[3]	# 0/0 + 0/1 + ./.
		pair_gts = (mat_gt, pat_gt)
		if pair_gts == (2, 2):
			# 両親合っている
			# 両親の尤度は(1-p)^2
			return self.log_likelihood_one(n+2, m)
		elif mat_gt == 2 or pat_gt == 2:
			# 片親だけ合っている
			# 両親の尤度は(1-p)p/3
			return self.log_likelihood_one(n+1, m+1)
		else:
			# 両親とも合っていない
			# 両親の尤度は(p/3)^2
			return self.log_likelihood_one(n, m+2)
	
	def determine(self, mat_gt: int, pat_gt: int,
						counter: list[int]) -> list[tuple[ParentComb, float]]:
		ls = [self.log_likelihood00(mat_gt, pat_gt, counter),
			  self.log_likelihood01(mat_gt, pat_gt, counter),
			  self.log_likelihood11(mat_gt, pat_gt, counter),
			  self.log_likelihood02(mat_gt, pat_gt, counter),
			  self.log_likelihood12(mat_gt, pat_gt, counter),
			  self.log_likelihood22(mat_gt, pat_gt, counter)]
		max_index, max_l = max(enumerate(ls), key=lambda p: p[1])
		pairs: list[tuple[ParentComb, float]] = [(ParentComb(max_index), max_l)]
		for i in range(6):
			if i != max_index and ls[i] > max_l - log(30):
				pairs.append((ParentComb(i), ls[i]))
		return pairs
