from __future__ import annotations

# coding: utf-8
# ImputeByRef.py
# リファレンスでimputeする

from itertools import *
import sys
from typing import Optional

from VCF import *
from GenoRecord import GenoRecord
from VCFGeno import VCFGeno
from pedigree import PedigreeTable, Family
import LargeFamilyRef
import LargeSelfFamilyRef
import SmallFamilyRef
from materials import *
from SampleManager import SampleManager
from Map import *
from option import *
from OptionSmall import OptionSmall
from common import *


#################### ImputeByRef ####################

def remove_reference_samples(vcf: VCFGeno, ref_samples: list[str]) -> VCFGeno:
	set_samples = set(ref_samples)
	samples = list(set(vcf.samples) - set_samples)
	cs: list[int] = [ c for c, s in enumerate(vcf.samples)
											if s not in set_samples ]
	records: list[GenoRecord] = []
	for record in vcf.records:
		geno = [ record.geno[c] for c in cs ]
		records.append(GenoRecord(record.pos, geno))
	return VCFGeno(samples, records, vcf.vcf)

def impute_vcf_chr(orig_vcf: VCFSmall, ref_vcf_: VCFSmall,
								sample_man: SampleManager,
								gmap: Map, option: Option) -> Optional[VCFGeno]:
	print('chr: %s %d records' % (orig_vcf.records[0].chrom(), len(orig_vcf)))
	sys.stdout.flush()
	
	ref_vcf = VCFGeno.convert(ref_vcf_)
	
	# 大きい家系
	merged_vcf = LargeFamilyRef.impute(sample_man.large_families, orig_vcf,
														ref_vcf, gmap, option)
	if merged_vcf is not None:
		sample_man.set(merged_vcf.samples)
	
	# 大きな自殖家系
	self_families = sample_man.extract_self_parent_non_imputed_families()
	imputed_vcf = LargeSelfFamilyRef.impute(self_families, orig_vcf,
											merged_vcf, ref_vcf, gmap, option)
	if imputed_vcf:
		merged_vcf = imputed_vcf
		sample_man.set(imputed_vcf.samples)
	
	# 小さい家系
	op_small = OptionSmall(gmap, option.num_threads, option.precision_ratio,
											option.imputes_isolated_samples,
											option.outputs_unimputed_samples)
	
	if merged_vcf is not None:
		sample_man.set(merged_vcf.samples)
	merged_vcf = SmallFamilyRef.impute(orig_vcf, merged_vcf,
											ref_vcf, op_small, sample_man,
											option.imputes_isolated_samples)
	merged_vcf.vcf = ref_vcf_	# writeのためにこれが必要
	
	sample_man.clear()
	return remove_reference_samples(merged_vcf, ref_vcf.samples)

def impute(vcf: VCFHuge, materials: Materials, option: Option) -> None:
	samples = vcf.samples
	ped = materials.ped.limit_samples(samples)
	ref_vcf = VCFHuge.read(option.path_ref_VCF)
	sample_man = SampleManager.create(ped, samples, ref_vcf.samples,
										option.lower_progs, option.families)
	sample_man.display_info(sys.stderr)
	iter = option.chroms_efficients()
	iter_vcf = vcf.divide_into_chromosomes()
	vcf_chr = next(iter_vcf)
	first = True
	# refには全部のchromがあるとする
	for b, ref_chr, gmap in zip(iter, ref_vcf.divide_into_chromosomes(),
												materials.iter_chr_maps()):
		# refにあってもvcfにないchromは飛ばす
		if vcf_chr.records[0].v[0] != ref_chr.records[0].v[0] or not b:
			try:
				vcf_chr = next(iter_vcf)
			except StopIteration:
				break
			continue
		
		vcf_imputed = impute_vcf_chr(vcf_chr, ref_chr,
											sample_man, gmap, option)
		if vcf_imputed is None:
			continue	# ここには来ないはず
		with open(option.path_out, 'w' if first else 'a') as out:
			vcf_imputed.write(out, with_header=first)
		first = False
		try:
			vcf_chr = next(iter_vcf)
		except StopIteration:
			break

__all__ = ['impute']
