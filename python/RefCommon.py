from __future__ import annotations

# coding: utf-8
# RefCommon.py

from VCF import VCFSmall
from GenoRecord import GenoRecord
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import VCFFamilyRecord
from KnownFamily import KnownFamily
from Genotype import Genotype
from typing import Optional


#################### RefCommon ####################

# phasedのrecordとnon-phasedのrecordをmergeする
# phasedにしかないpositionは./.で埋めて
# phasedにないpositionは捨てる
def merge_records_core(phased_vcf: VCFGenoBase,
						non_phased_vcf: VCFGenoBase,
						samples: list[str]) -> list[tuple[int, list[int]]]:
	# 全部N/Aで埋める
	N = len(samples)
	M = len(phased_vcf)
	poss = [ phased_vcf.get_record(i).pos for i in range(M) ]
	genos: list[tuple[int, list[int]]] = [ (pos, [Genotype.NA]*N)
														for pos in poss ]
	cs1 = phased_vcf.extract_columns(samples)
	cs2 = non_phased_vcf.extract_columns(samples)
	
	# phasedのサンプルのGenotypeを入れる
	for j in range(M):
		record1 = phased_vcf.get_record(j)
		for i, c in enumerate(cs1):
			if c != -1:
				genos[j][1][i] = record1.geno[c]
	
	k, l = 0, 0
	L1, L2 = len(phased_vcf), len(non_phased_vcf)
	while k < L1 and l < L2:
		record1 = phased_vcf.get_record(k)
		record2 = non_phased_vcf.get_record(l)
		if record1.pos == record2.pos:
			for i, c in enumerate(cs2):
				if c != -1 and cs1[i] == -1:
					genos[k][1][i] = record2.geno[c]
			k += 1
			l += 1
		elif record1.pos < record2.pos:
			k += 1
		else:
			l += 1
	
	return genos

def merge_records(phased_vcf: VCFGenoBase,
						non_phased_vcf: VCFGenoBase,
						samples: list[str]) -> list[GenoRecord]:
	genos = merge_records_core(phased_vcf, non_phased_vcf, samples)
	return [ GenoRecord(pos, geno) for pos, geno in genos ]

def merge_family_records(phased_vcf: VCFGenoBase,
						non_phased_vcf: VCFGenoBase,
						samples: list[str]) -> list[VCFFamilyRecord]:
	genos = merge_records_core(phased_vcf, non_phased_vcf, samples)
	return [ VCFFamilyRecord(pos, geno) for pos, geno in genos ]

# phased_vcfのポジションでvcfに無いところはN/Aにする
def expand_records(vcf: VCFGeno, phased_vcf: VCFGeno) -> list[GenoRecord]:
	if len(vcf) == 0:
		return []
	
	N = len(vcf.samples)
	M = len(phased_vcf)
	poss = [ phased_vcf.records[i].pos for i in range(M) ]
	records = [ GenoRecord(pos, [Genotype.NA]*N) for pos in poss ]
	
	k, l = 0, 0
	L1, L2 = len(vcf), len(records)
	while k < L1 and l < L2:
		record1 = vcf.records[k]
		record2 = records[l]
		if record1.pos == record2.pos:
			for i in range(N):
				record2.geno[i] = record1.geno[i]
			k += 1
			l += 1
		elif record1.pos < record2.pos:
			k += 1
		else:
			l += 1
	
	return records
