# coding: utf-8
# impute_prog_only.py

from __future__ import annotations
from itertools import *
import sys
from typing import Optional

from VCF import *
from VCFGeno import VCFGenoBase
from VCFFamily import VCFFamily, VCFFamilyRecord
from VCFImpFamilyRecord import FillType, VCFImpFamilyRecord
from VCFHeteroHomo import VCFHeteroHomoRecord
from VCFHeteroHomoPP import *
from VCFHomoHomo import VCFHomoHomoRecord
from VCFHeteroHeteroLiteRecord import VCFHeteroHeteroLiteRecord
from VCFJunkRecord import VCFJunkRecord
from VCFFillableRecord import VCFFillableRecord
from VCFFillable import VCFFillable
import SmallFamily
from pedigree import PedigreeTable, Family
from SampleManager import SampleManager
import ClassifyRecord as CR
from TypeDeterminer import ParentComb, TypeDeterminer
from Map import *
from option import *
from common import *


#################### process ####################

def classify_record(i: int, record: VCFFamilyRecord,
					samples: list[str], td: TypeDeterminer,
					heho_records: list[Optional[VCFHeteroHomoRecord]],
					other_records: list[Optional[VCFImpFamilyRecord]]) -> None:
	wrong_type, pair = CR.classify_record(record, td, False)
	pos = record.pos
	geno = record.geno
	if pair.is_homohomo():
		record_ = VCFHomoHomoRecord(pos, geno, i, wrong_type, pair)
		other_records[i] = record_.impute()
	elif pair.is_heterohomo():
		heho_records[i] = VCFHeteroHomoRecord(pos, geno, i, wrong_type, pair)
	elif pair == ParentComb.P01x01:
		other_records[i] = VCFHeteroHeteroLiteRecord(pos, geno, i, wrong_type)
	else:		# 候補が無い
		other_records[i] = VCFJunkRecord(pos, geno, i, wrong_type)

def classify_records(vcf: VCFFamily, option: Option
									) -> tuple[list[VCFHeteroHomoRecord],
												 list[VCFImpFamilyRecord]]:
	td = CR.get_typedeterminer(vcf.num_samples()-2, option.ratio)
	# 整理したい
	records1: list[Optional[VCFHeteroHomoRecord]] = [None] * len(vcf)
	records2: list[Optional[VCFImpFamilyRecord]] = [None] * len(vcf)
	for i, record in enumerate(vcf.records):
		classify_record(i, record, vcf.samples, td, records1, records2)
	
	heho_records: list[VCFHeteroHomoRecord] = []
	other_records: list[VCFImpFamilyRecord] = []
	for i in range(len(vcf)):
		r1 = records1[i]
		r2 = records2[i]
		if r1 is not None:
			heho_records.append(r1)
		elif r2 is not None:
			other_records.append(r2)
	
	return (heho_records, other_records)

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
	# GTのみにする
	header = vcf_parents.trim_header(samples)
	records: list[VCFRecord] = []
	j = 0
	for i, record1 in enumerate(vcf_parents.records):
		if j == len(vcf_progenies):
			record = fill_NA(record1, samples)
		else:
			record2 = vcf_progenies.records[j]
			if record1.pos == record2.pos:
				record = merge_record(record1, record2, samples)
				j += 1
			else:
				record = fill_NA(record1, samples);
		records.append(record)
	return VCFSmall(header, records)

def merge_vcf(rss: list[list[VCFFillableRecord]], samples: list[str],
								gmap: Map, vcf: VCFSmall) -> VCFHeteroHomoPP:
	rs = rss[0] + rss[1] + rss[2] + rss[3]
	rs.sort(key=lambda r: r.pos)
	return VCFHeteroHomoPP(samples, rs, gmap, vcf)

def impute_prog_vcf_chr(parent_vcf: VCFSmall, prog_vcf: VCFSmall,
									gmap: Map, option: Option) -> VCFGenoBase:
	print('chr: %s %d records' % (parent_vcf.records[0].chrom(),
													len(parent_vcf)))
	sys.stdout.flush()
	# parent_vcfは両親のみ、prog_vcfは後代のみという前提
	samples = parent_vcf.samples + prog_vcf.samples
	orig_vcf = merge_parents_progenies(parent_vcf, prog_vcf, samples)
	merged_vcf = VCFHeteroHomoPP.merge(parent_vcf, prog_vcf, orig_vcf,
														samples, gmap, option)
	rss = VCFHeteroHomoPP.classify_records(samples, merged_vcf.records,
																orig_vcf)
	mat_vcf = VCFHeteroHomoPP(samples, rss[FillType.MAT.value], gmap, orig_vcf)
	pat_vcf = VCFHeteroHomoPP(samples, rss[FillType.PAT.value], gmap, orig_vcf)
	mat_vcf.impute()
	pat_vcf.impute()
	merged_vcf = merge_vcf(rss, samples, gmap, orig_vcf)
	merged_vcf.fill()
	return merged_vcf

__all__ = ['impute_prog_vcf_chr']
