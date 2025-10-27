from __future__ import annotations

# coding: utf-8
# LargeSelfFamily.py

from itertools import *
from collections import defaultdict
import sys

from typing import Optional

from VCF import *
from pedigree import PedigreeTable, Family
from KnownFamily import KnownFamily
from VCFGeno import VCFGeno
from GenoRecord import GenoRecord
from VCFSelfHetero import *
from VCFSelfHeteroRecord import VCFSelfHeteroRecord
from VCFSelfHomoRecord import *
from VCFSelfParentImputed import VCFSelfParentImputed
from VCFSelfJunkRecord import VCFSelfJunkRecord
from VCFSelfFillable import *
from VCFImpSelfRecord import VCFImpSelfRecord
from SampleManager import SampleManager
from Map import *
import ClassifyRecord as CR
from TypeDeterminer import TypeDeterminer, ParentComb
from option import *
from common import *


#################### process ####################

def divide_records(vcf: VCFGeno, op: Option
								) ->tuple[list[VCFSelfHeteroRecord],
										  list[VCFImpSelfRecord]]:
	td = CR.get_typedeterminer(vcf.num_samples()-1, op.ratio)
	he_records: list[VCFSelfHeteroRecord] = []
	other_records: list[VCFImpSelfRecord] = []
	samples = vcf.samples
	for i, record in enumerate(vcf.records):
		pos = record.pos
		geno = record.geno
		wrong_type, pair = CR.classify_self_record(record, td)
		if pair == ParentComb.P01x01:
			record1 = VCFSelfHeteroRecord(pos, geno, i, wrong_type, pair)
			he_records.append(record1)
		elif pair.is_homohomo():
			record2 = VCFSelfHomoRecord(pos, geno, i, wrong_type, pair)
			other_records.append(record2.impute())
		else:
			other_records.append(VCFSelfJunkRecord(pos, geno, i, wrong_type))
	
	return (he_records, other_records)

def extract_parents(vcfs: list[VCFSelfFillable]) -> VCFGeno:
	samples = [ vcf.get_samples()[0] for vcf in vcfs ]
	M = len(vcfs[0])
	records: list[GenoRecord] = []
	for i in range(M):
		pos = vcfs[0].get_record(i).pos
		geno = [ vcf.get_record(i).geno[0] for vcf in vcfs ]
		record = GenoRecord(pos, geno)
		records.append(record)
	return VCFGeno(samples, records, vcfs[0].vcf)

def impute(families: list[KnownFamily], orig_vcf: VCFSmall,
		   merged_vcf: VCFGeno, geno_map: Map, op: Option) -> Optional[VCFGeno]:
	if not families:
		return None
	
	vcfs: list[VCFSelfFillable] = []
	for family in families:
		samples = [family.mat] + family.progenies
		vcf = VCFGeno.extract_samples(samples, orig_vcf)
		he_records, other_records = divide_records(vcf, op)
		vcf_hetero = VCFSelfHetero(samples, he_records, geno_map, orig_vcf)
		vcf_heteros, unused = vcf_hetero.impute()
		other_records.extend(unused)
		vcf_filled = VCFSelfFillable.fill(vcf_heteros, other_records)
		vcfs.append(vcf_filled)
	
	vcf_parents = extract_parents(vcfs)
	merged_vcf = VCFGeno.join([merged_vcf, vcf_parents], orig_vcf.samples)
	return merged_vcf

__all__ = ['impute_large_self_families']
