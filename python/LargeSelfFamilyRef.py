from __future__ import annotations

# coding: utf-8
# LargeSelfFamilyRef.py

from itertools import *
from collections import defaultdict
import sys

from typing import Optional

from VCF import *
from pedigree import PedigreeTable, Family
from KnownFamily import KnownFamily
from GenoRecord import GenoRecord
from VCFGeno import VCFGenoBase, VCFGeno
from VCFSelfImputable import *
from VCFSelfParentImputed import VCFSelfParentImputed
from VCFSelfProgenyImputed import VCFSelfProgenyImputed
from VCFSelfImputable import *
import RefCommon
from Map import *
from option import *
from Genotype import Genotype
from common import *


#################### LargeSelfFamilyRef ####################

def create_parent_phased_vcf(family: KnownFamily, orig_vcf: VCFSmall,
								phased_vcf: VCFGeno, gmap: Map,
								op: Option) -> Optional[VCFSelfParentImputed]:
	# phasedの後代を取り除いたsamplesを作る
	ids = phased_vcf.extract_columns(family.samples())	# phased_vcf内のindex
	prog_ids = [ i for i in range(1, len(family.samples())) if ids[i] == -1 ]
	if len(prog_ids) == 0:
		return None
	
	samples = [family.mat] + [ family.samples()[i] for i in prog_ids ]
	vcf_progs = VCFGeno.extract_samples(samples[1:], orig_vcf)
	records = RefCommon.merge_records(phased_vcf, vcf_progs, samples)
	return VCFSelfParentImputed(samples, records, gmap, phased_vcf.vcf)

def create_progeny_phased_vcf(family: KnownFamily, orig_vcf: VCFSmall,
								ref_vcf: VCFGeno, phased_vcf: VCFGeno,
								gmap: Map, op: Option) -> VCFSelfProgenyImputed:
	# phasedの後代を一つだけ残したsamplesを作る
	N = len(family.samples())
	ids = phased_vcf.extract_columns(family.samples())	# phased_vcf内のindex
	prog_ids = [ i for i in range(1, len(family.samples())) if ids[i] != -1 ]
	phased_prog_id = prog_ids[0]
	phased_sample = family.samples()[phased_prog_id]
	set_prog_ids = set(prog_ids[1:])	# 最初の1個だけ残す
	samples = [ family.samples()[i] for i in range(N) if i not in set_prog_ids ]
	non_phased_samples = [ s for s in samples if s != phased_sample ]
	
	non_phased_vcf = VCFGeno.extract_samples(non_phased_samples, orig_vcf)
	records = RefCommon.merge_records(phased_vcf, non_phased_vcf, samples)
	ref_haps = ref_vcf.create_ref_haps()
	return VCFSelfProgenyImputed(samples, records, ref_haps,
										prog_ids[0], gmap, phased_vcf.vcf)

def impute(families: list[KnownFamily], orig_vcf: VCFSmall,
								merged_vcf: Optional[VCFGeno],
								ref_vcf: VCFGeno,
								gmap: Map, op: Option) -> Optional[VCFGeno]:
	if not families:
		return None
	
	# 今までimputeしたVCFとリファレンスを合わせてphasedという
	phased_vcf = VCFGeno.merge(merged_vcf, ref_vcf)
	
	vcfs: list[VCFGenoBase] = []
	set_phased = set(phased_vcf.samples)
	for family in families:
		if family.mat in set_phased:
			# 親がphasedならば、そのままimputeする
			vcf1 = create_parent_phased_vcf(family, orig_vcf,
												phased_vcf, gmap, op)
			if vcf1 is not None:
				vcf1.impute()
				vcfs.append(vcf1)
		elif set_phased & set(family.samples()):
			# 後代にphasedがあれば、そのままimputeする
			vcf2 = create_progeny_phased_vcf(family, orig_vcf, ref_vcf,
														phased_vcf, gmap, op)
			vcf2.impute()
			vcfs.append(vcf2)
	
	if merged_vcf is not None:
		vcfs.append(merged_vcf)
	
	if vcfs:
		return VCFGeno.join(vcfs, orig_vcf.samples)
	else:
		return None

__all__ = ['impute']
