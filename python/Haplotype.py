from __future__ import annotations

# coding: utf-8
# VCFImputable.py

from itertools import *
from abc import ABC, abstractmethod
from typing import Iterator, Tuple


#################### Haplotype ####################

class Haplotype(object):
	def __init__(self, h: list[int], sample_id: int, i: int):
		self.hap: list[int] = h
		self.position: tuple[int, int] = (sample_id, i)
	
	@staticmethod
	def default_value() -> Haplotype:
		hap: list[int] = []
		return Haplotype(hap, 0, 2)
	
	@staticmethod
	def score(int_gts: list[int], hap_mat: Haplotype,
								  hap_pat: Haplotype) -> int:
		return sum(1 for h1, h2, gt in zip(hap_mat.hap, hap_pat.hap, int_gts)
															if h1 + h2 == gt)
	
	@staticmethod
	def collect_optimal_haplotype_pairs(int_gts: list[int],
										haps_mat: list[Haplotype],
										haps_pat: list[Haplotype]
												) -> list[HaplotypePair]:
		scores = [ Haplotype.score(int_gts, hap_mat, hap_pat)
						for hap_mat, hap_pat in product(haps_mat, haps_pat) ]
		max_score = max(scores)
		return [ hap_pair
					for s, hap_pair in zip(scores, product(haps_mat, haps_pat))
															if s == max_score ]
	
	@staticmethod
	def match_score(prev_hap: HaplotypePair, hap: HaplotypePair) -> int:
		prev_hap1, prev_hap2 = prev_hap
		hap1, hap2 = hap
		return ((1 if prev_hap1.position == hap1.position else 0) +
				(1 if prev_hap2.position == hap2.position else 0))
	
	@staticmethod
	def collect_max_score(combs: list[HaplotypePair],
								prev_hap: HaplotypePair) -> list[HaplotypePair]:
		scores = [ (Haplotype.match_score(prev_hap, hap), hap)
													for hap in combs ]
		max_score = max(s for s, hap in scores)
		return [ hap for s, hap in scores if s == max_score ]
	
	@staticmethod
	def impute(int_gts: list[int], haps_mat: list[Haplotype],
						haps_pat: list[Haplotype], prev_hap: HaplotypePair,
						seed: int) -> HaplotypePair:
		# blute-force
		combs = Haplotype.collect_optimal_haplotype_pairs(int_gts,
														haps_mat, haps_pat)
		
		# collect combinations that are most matched with the previous
		filtered_combs = Haplotype.collect_max_score(combs, prev_hap)
		
		# select a combination with seed
		j = seed % len(filtered_combs)
		return filtered_combs[j]
	
HaplotypePair = Tuple[Haplotype, Haplotype]
