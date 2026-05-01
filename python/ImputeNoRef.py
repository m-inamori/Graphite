from __future__ import annotations

# coding: utf-8
# ImputeNoRef.py
# リファレンス無しでimputeする

from itertools import *
import sys
from typing import Iterator, Optional

from VCF import *
from VCFGeno import VCFGeno
from pedigree import PedigreeTable, Family
import LargeFamily
import LargeSelfFamily
import SmallFamily
from materials import *
from SampleManager import SampleManager
from Map import *
from option import *
from OptionSmall import OptionSmall
from common import *


#################### ImputeNoRef ####################

def impute_vcf_chr(orig_vcf: VCFSmall, sample_man: SampleManager,
							geno_map: Map, option: Option) -> Optional[VCFGeno]:
	print('chr: %s %d records' % (orig_vcf.records[0].chrom(), len(orig_vcf)))
	sys.stdout.flush()
	merged_vcf = LargeFamily.impute(orig_vcf, sample_man.large_families,
															geno_map, option)
	if merged_vcf is not None:
		sample_man.set(merged_vcf.samples)
	
	self_families = sample_man.extract_self_parent_non_imputed_families()
	imputed_vcf = LargeSelfFamily.impute(self_families, orig_vcf, merged_vcf,
															geno_map, option)
	if imputed_vcf:
		merged_vcf = imputed_vcf
		sample_man.set(imputed_vcf.samples)
	
	op_small = OptionSmall(geno_map, option.num_threads, option.precision_ratio,
											option.imputes_isolated_samples,
											option.outputs_unimputed_samples)
	
	# 補完できる家系がなくなるまで繰り返す
	if merged_vcf:
		sample_man.set(merged_vcf.samples)
		merged_vcf = SmallFamily.impute(orig_vcf, merged_vcf,
												op_small, sample_man,
												option.imputes_isolated_samples)
	
	sample_man.clear()
	return merged_vcf

def impute(vcf: VCFHuge, materials: Materials, option: Option) -> None:
	samples = vcf.samples
	ped = materials.ped.limit_samples(samples)
	sample_man = SampleManager.create(ped, samples, [],
										option.lower_progs, option.families)
	sample_man.display_info(sys.stderr)
	
	iter = option.chroms_efficients()
	first = True
	for b, vcf_chr, gmap in zip(iter, vcf.divide_into_chromosomes(),
													materials.iter_chr_maps()):
		if not b:
			continue
		vcf_chr.check_records()
		vcf_imputed = impute_vcf_chr(vcf_chr, sample_man, gmap, option)
		if vcf_imputed is None:
			continue
		with open(option.path_out, 'w' if first else 'a') as out:
			vcf_imputed.write(out, with_header=first)
		first = False

__all__ = ['impute']
