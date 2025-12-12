from __future__ import annotations

# coding: utf-8
# ReferenceHaplotype.py

from itertools import *
from collections import Counter

from VCFGeno import VCFGenoBase, VCFGeno
from SampleManager import *
from Genotype import Genotype
from common import *


#################### ReferenceHaplotype ####################

def extract_haplotypes(phased_vcf: VCFGeno,
						sample_man: SampleManager) -> list[list[int]]:
	reference: list[str] = sample_man.collect_reference()
	ref_columns = phased_vcf.extract_columns(reference)
	
	N = len(reference)
	M = len(phased_vcf)
	
	hs = [ h for c in ref_columns for h in range(c*2, c*2+2) ]
	gts = [ [ record.get_allele(h>>1, h&1) for record in phased_vcf.records ]
																for h in hs ]
	MIN_REF_NUM = 10
	if len(gts) <= MIN_REF_NUM:
		return gts
	
	K = min(10, M)
	# ローリングハッシュ
	a: list[list[int]] = [ [0] * (M-K+1) for _ in range(N*2) ]
	for h in range(N*2):
		n = 0
		for i in range(M):
			gt = gts[h][i]
			if i < K - 1:
				n = (n << 1) | gt
			else:
				n = ((n << 1) & ((1 << 10) - 1)) | gt
				a[h][i-K+1] = n
	
	# 同じハッシュ値があったらどのハプロタイプと同じか調べる
	b: list[list[int]] = [ [0] * (M-K+1) for _ in range(N*2) ]
	for i in range(M-K+1):
		for h in range(N*2):
			for l in range(h):
				if a[l][i] == a[h][i]:
					b[h][i] = l
					break
			else:
				b[h][i] = h
	
	# Genotypeに戻ってできるだけ自分以外になるようにする
	for k in range(1, N*2):
		i = 0
		c: list[tuple[int, int, int]] = []	# [(index, first, last)]
		for h, v1 in groupby(b[k]):
			next_i = i + sum(1 for _ in v1)
			c.append((h, i, next_i))
			i = next_i
		
		# 自分以外の前後を可能なら拡張する
		for j in range(len(c)):
			h0, first0, last0 = c[j]
			if h0 != k:
				continue
			if j != 0:
				# 前を見る
				h1, first1, last1 = c[j-1]
				for i in range(first0, last0):
					if gts[h1][i] != gts[k][i]:
						break
					b[k][i] = h1
			if j != len(c) - 1:
				# 後ろを見る
				h2, first2, last2 = c[j+1]
				for i in range(last0-1, first0, -1):
					if gts[h2][i] != gts[k][i]:
						break
					b[k][i] = h2
	
	# 2つのハプロタイプが長さの1割以下しかないサンプルは捨てる
	# ただし、MIN_REF_NUM以下にならないようにする
	counter = Counter(g for v in b for g in v)
	w = sorted((counter[i], i) for i in range(N*2))
	if w[N*2-MIN_REF_NUM][0] * 10 >= M:
		return [ gts[i] for c, i in w if c * 10 >= M ]
	else:
		return [ gts[i] for c, i in w[N*2-MIN_REF_NUM:] ]

def count_same_alleles(gts: list[int], ref_hap: list[int]) -> int:
	n = 0
	for k in range(len(gts)):
		gt = Genotype.unphased(gts[k])
		if gt == Genotype.UN_00:
			if ref_hap[k] == 0:
				n += 1
#		elif gt == Genotype.UN_01:
#			n += 1
		elif gt == Genotype.UN_11:
			if ref_hap[k] == 1:
				n += 1
	return n

# lowerを下回らないようにgtsに似た
def filter_haplotypes(ref_haps: list[list[int]],
						gts: list[int], lower: int) -> list[list[int]]:
	# 単純に染色体全部のホモを見て、同じアレルの数を数える
	NH = len(ref_haps)
	# [(同じアレルの数, ref_hapsのindex)]
	ns = [ (count_same_alleles(gts, ref_haps[i]), i) for i in range(NH) ]
	ns.sort()
	return [ ref_haps[i] for n, i in ns[-lower:] ]

__all__ = ['extract_haplotypes', 'filter_haplotypes']
