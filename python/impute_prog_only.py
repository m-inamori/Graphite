# coding: utf-8
# impute_prog_only.py

from __future__ import annotations
from itertools import *
import sys
from typing import Optional

from VCF import *
from VCFGeno import VCFGenoBase, VCFGeno
from VCFBothParentImputed import VCFBothParentImputed
from VCFFamily import VCFFamily, VCFFamilyRecord
import SmallFamily
from pedigree import PedigreeTable, Family
from Map import *
from option import *
from common import *


#################### process ####################

def fill_NA(record1: VCFRecord, samples: list[str]) -> VCFRecord:
	NA_len = len(samples) - len(record1.samples)
	parents_v = record1.v[:8] + ['GT'] + [ gt[:3] for gt in record1.v[9:] ]
	v = parents_v + ['./.'] * NA_len
	return VCFRecord(v, samples)

def merge_record(record1: VCFRecord, record2: VCFRecord,
										samples: list[str]) -> VCFRecord:
	v = record1.v + record2.v[9:]
	v[8] = 'GT'
	for c in range(9, len(v)):
		v[c] = v[c][:3]
	return VCFRecord(v, samples)

def merge_parents_progenies(vcf_parents: VCFSmall, vcf_progenies: VCFSmall,
												samples: list[str]) -> VCFSmall:
	# 後代で無いポジションはN/Aで埋める
	header = vcf_parents.trim_header(samples)
	records: list[VCFRecord] = []
	j = 0
	for i, record1 in enumerate(vcf_parents.records):
		if j == len(vcf_progenies):
			record = fill_NA(record1, samples)
		else:
			record2 = vcf_progenies.records[j]
			if record1.pos() == record2.pos():
				record = merge_record(record1, record2, samples)
				j += 1
			else:
				record = fill_NA(record1, samples);
		records.append(record)
	return VCFSmall(header, records)

def impute_prog_vcf_chr(parent_vcf: VCFSmall, prog_vcf: VCFSmall,
									gmap: Map, option: Option) -> VCFGenoBase:
	print('chr: %s %d records' % (parent_vcf.records[0].chrom(),
													len(parent_vcf)))
	sys.stdout.flush()
	# parent_vcfは両親のみ、prog_vcfは後代のみという前提
	samples = parent_vcf.samples + prog_vcf.samples
	orig_vcf = merge_parents_progenies(parent_vcf, prog_vcf, samples)
	merged_vcf = VCFGeno.convert(orig_vcf)
	records = [ VCFFamilyRecord(r.pos, r.geno) for r in merged_vcf.records ]
	vcf = VCFBothParentImputed(merged_vcf.samples, records, gmap, orig_vcf)
	vcf.impute(1.0, 1)
	imputed_prog_vcf = vcf.extract_by_samples(prog_vcf.samples)
	return imputed_prog_vcf

__all__ = ['impute_prog_vcf_chr']
