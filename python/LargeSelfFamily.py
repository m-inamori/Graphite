from __future__ import annotations

# coding: utf-8
# LargeFamily.py

from itertools import *
from collections import defaultdict
import sys

from typing import Optional

from VCF import *
from pedigree import PedigreeTable, Family
from KnownFamily import KnownFamily
from VCFSelfHetero import *
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

def divide_records(vcf: VCFSmall, op: Option
								) ->tuple[list[VCFSelfHeteroRecord],
										  list[VCFImpSelfRecord]]:
	td = CR.get_typedeterminer(vcf.num_samples()-1, op.ratio)
	he_records: list[VCFSelfHeteroRecord] = []
	other_records: list[VCFImpSelfRecord] = []
	samples = vcf.samples
	for i, record in enumerate(vcf.records):
		v = record.v
		wrong_type, pair = CR.classify_self_record(record, td)
		if pair == ParentComb.P01x01:
			record1 = VCFSelfHeteroRecord(v, samples, i, wrong_type, pair)
			he_records.append(record1)
		elif pair.is_homohomo():
			record2 = VCFSelfHomoRecord(v, samples, i, wrong_type, pair)
			other_records.append(record2.impute())
		else:
			other_records.append(VCFSelfJunkRecord(v, samples, i, wrong_type))
	
	return (he_records, other_records)

def extract_parents(vcfs: list[VCFSelfFillable]) -> VCFSmall:
	samples = [ vcf.get_samples()[0] for vcf in vcfs ]
	header = vcfs[0].trim_header(samples)
	M = len(vcfs[0])
	records: list[VCFRecord] = []
	for i in range(M):
		v = vcfs[0].get_record(i).v[:10]
		for j in range(1, len(vcfs)):
			v.append(vcfs[j].get_record(i).v[9])
		record = VCFRecord(v, samples)
		records.append(record)
	return VCFSmall(header, records)

def impute(orig_vcf: VCFSmall, merged_vcf: VCFSmall,
				sample_man: SampleManager, geno_map: Map, op: Option
									) -> Optional[VCFSmallBase]:
	families = sample_man.extract_self_parent_non_imputed_families()
	if not families:
		return None
	
	vcfs = []
	for family in families:
		samples = [family.mat] + family.progenies
		vcf = orig_vcf.extract_samples(samples)
		he_records, other_records = divide_records(vcf, op)
		vcf_hetero = VCFSelfHetero(vcf.header, he_records, geno_map)
		vcf_heteros, unused = vcf_hetero.impute()
		other_records.extend(unused)
		vcf_filled = VCFSelfFillable.fill(vcf_heteros, other_records)
		vcfs.append(vcf_filled)
	
	vcf_parents = extract_parents(vcfs)
	merged_vcf = VCFSmall.join([merged_vcf, vcf_parents], orig_vcf.samples)
	return merged_vcf

__all__ = ['impute_large_self_families']
