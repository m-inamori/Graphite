from __future__ import annotations

# coding: utf-8
# Imputer.py

from itertools import *
from collections import defaultdict
from math import exp

from typing import Generator, Dict

import Baum_Welch_with_fixed_Ts_log as BWTL
from common import *


#################### Region ####################

class Region:
	def __init__(self, color: str, first: int, last: int, cM: float):
		self.color: str	= color
		self.first: int	= first
		self.last: int	= last
		self.cM: float	= cM
	
	@staticmethod
	def create(seq: str, cMs: list[float]) -> list[Region]:
		regions = []
		last = 0
		for color, w in groupby(seq):
			n = sum(1 for _ in w)
			first, last = last, last + n
			cM1 = cMs[0] if first == 0 else (cMs[first-1] + cMs[first]) / 2
			cM2 = cMs[-1] if last == len(seq) else (cMs[last-1] + cMs[last]) / 2
			regions.append(Region(color, first, last, cM2 - cM1))
		return regions


#################### State ####################

class State:
	def __init__(self, regions: list[Region], colors, size: int,
						n: int, length: float, p_len: float, min_c: float):
		self.regions: list[Region] = regions
		# 効率化のためのリスト構造
		# typingでうまく型を表せない
		self.new_colors = colors	# [[[...], str], str]
		self.size: int = size
		self.num_continuous: int = n
		self.continuous_length: float = length
		self.painted_length: float = p_len
		self.MIN_CROSSOVER: float = min_c
	
	def last_color(self) -> str:
		return self.new_colors[1]
	
	def to_list(self) -> list[str]:
		cs = self.new_colors
		v: list[str] = []
		while cs:
			v = cs[1:] + v
			cs = cs[0]
		return v
	
	def extend(self, cM: float) -> State:
		painted = self.last_color() != self.next_color()
		return State(self.regions,
					 [self.new_colors, self.last_color()],
					 self.size + 1,
					 self.num_continuous + 1,
					 self.continuous_length + cM,
					 self.painted_length + (cM if painted else 0.0),
					 self.MIN_CROSSOVER)
	
	def change(self, cM: float) -> State:
		color = '1' if self.last_color() == '0' else '0'
		painted = color != self.next_color()
		return State(self.regions,
					 [self.new_colors, color],
					 self.size + 1,
					 1,
					 cM,
					 self.painted_length + (cM if painted else 0.0),
					 self.MIN_CROSSOVER)
	
	def current_color(self) -> str:
		return self.regions[self.size-1].color
	
	def next_color(self) -> str:
		return self.regions[self.size].color
	
	def is_extendable(self) -> bool:
		return (self.size != len(self.regions) - 1 or
				self.next_color() == self.last_color())
	
	def rev_colors(self) -> Generator[str, None, None]:
		new_colors = self.new_colors
		while new_colors:
			yield new_colors[1]
			new_colors = new_colors[0]
	
	def is_changable(self) -> bool:
		color = self.next_color()
		return (self.continuous_length >= self.MIN_CROSSOVER or
				# it is painted the same color till now
				all(c == self.last_color() for c in self.rev_colors()))
	
	def nexts(self, cM: float) -> Generator[State, None, None]:
		if self.is_extendable():
			yield self.extend(cM)
		if self.is_changable():
			yield self.change(cM)
	
	@staticmethod
	def init(regions: list[Region], min_c: float) -> State:
		new_colors = [[], regions[0].color]
		return State(regions, new_colors, 1, 1, regions[0].cM, 0.0, min_c)


#################### paint ####################

def paint(seq: str, cMs: list[float], min_c: float) -> str:
	def select(states: list[State]) -> list[State]:
		dic: Dict[tuple[int, str], list[State]] = defaultdict(list)
		for state in states:
			dic[(state.num_continuous, state.last_color())].append(state)
		
		return [ min(ss, key=lambda s: s.painted_length)
									for ss in dic.values() ]
	
	regions = Region.create(seq, cMs)
	if len(regions) <= 2:
		return seq
	
	# DP
	states: list[State] = [State.init(regions, min_c)]
	for region in regions[1:]:
		new_states = [ s2 for s1 in states for s2 in s1.nexts(region.cM) ]
		states = select(new_states)
	
	state = min(states, key=lambda s: s.painted_length)
	new_colors = list(state.rev_colors())[::-1]
	return ''.join(color * (region.last - region.first)
						for region, color in zip(regions, new_colors))


#################### impute ####################

def is_all_same_without_N(seq: str) -> bool:
	prev_c = '-'
	for c in seq:
		if c == 'N':
			continue
		elif prev_c == '-':
			prev_c = c
		elif c != prev_c:
			return False
	else:
		return True

def create_same_color_string(seq: str, default_color: str):
	if all(c == 'N' for c in seq):
		c_base = default_color
	else:
		c_base = next(c for c in seq if c != 'N')
	return c_base * len(seq)

def compute_T(p: float, hidden_states: list[str]
							) -> Dict[tuple[str, str], float]:
	return { (h1, h2): 1.0-p if h1 == h2 else p
					for h1, h2 in product(hidden_states, repeat=2) }

def impute(seq: str, hidden_states: list[str],
					states: list[str], cMs: list[float]) -> str:
	ds = [ (p2-p1)/100 for p1, p2 in zip(cMs, cMs[1:]) ]
	ps = [ (exp(d) - exp(-d)) / (exp(d) + exp(-d)) for d in ds ]
	Ts = [ compute_T(p, hidden_states) for p in ps ]
	pi, E = BWTL.Baum_Welch(Ts, seq, states, hidden_states)
	A = BWTL.reverse_probs(E, hidden_states, states)
	hidden_seq = BWTL.Viterbi(seq, states, hidden_states, pi, Ts, A)
	return hidden_seq


#################### main ####################

__all__ = ['paint', 'impute', 'is_all_same_without_N',
								'create_same_color_string']
