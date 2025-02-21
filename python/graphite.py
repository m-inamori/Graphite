from __future__ import annotations

# coding: utf-8
# graphite.py
# フルでphasingする

from itertools import *
import sys
from typing import Iterator

from VCF import *
from pedigree import PedigreeTable, Family
from LargeFamily import correct_large_family_VCFs
import SmallFamily
from impute_prog_only import *
from materials import *
from SampleManager import SampleManager
from Map import *
from option import *
from common import *
from exception_with_code import *


#################### process ####################

def impute_vcf_chr(orig_vcf: VCFSmall, sample_man: SampleManager,
						geno_map: Map, option: Option) -> VCFSmall:
	print('chr: %s %d records' % (orig_vcf.records[0].chrom(), len(orig_vcf)))
	sys.stdout.flush()
	merged_vcf = correct_large_family_VCFs(orig_vcf, sample_man.large_families,
														geno_map, option)
	if option.only_large_families:
		return merged_vcf
	
	# 補完できる家系がなくなるまで繰り返す
	sample_man.set(merged_vcf.samples)
	
	merged_vcf = SmallFamily.impute_small_family_VCFs(orig_vcf, merged_vcf,
														geno_map, sample_man,
														option.num_threads)
	
	# 最後に孤立したサンプルを補完する
	samples = sample_man.extract_isolated_samples()
	if samples:
		new_imputed_vcf = SmallFamily.impute_isolated_samples(orig_vcf,
												merged_vcf, sample_man, samples,
												geno_map, option.num_threads)
		vcfs: list[VCFSmallBase] = [merged_vcf, new_imputed_vcf]
		merged_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	
	sample_man.clear()
	return merged_vcf

# 処理する染色体かどうかを排出する
def chroms_efficients(option: Option) -> Iterator[bool]:
	if not option.chroms:
		return repeat(True)
	
	upper = max(option.chroms)
	v: list[bool] = [False] * (upper + 1)
	for i in option.chroms:
		v[i] = True
	return (b for b in v)

def impute_all(vcf: VCFHuge, materials: Materials, option: Option) -> None:
	samples = vcf.samples
	ped = materials.ped.limit_samples(samples)
	sample_man = SampleManager.create(ped, samples,
										option.lower_progs, option.families)
	sample_man.display_info(sys.stderr)
	
	iter = chroms_efficients(option)
	first = True
	for b, vcf_chr, gmap in zip(iter, vcf.divide_into_chromosomes(),
													materials.iter_chr_maps()):
		if not b:
			continue
		vcf_chr.check_records()
		vcf_imputed = impute_vcf_chr(vcf_chr, sample_man, gmap, option)
		with open(option.path_out, 'w' if first else 'a') as out:
			vcf_imputed.write(out, with_header=first)
		first = False

def impute_progenies(vcf: VCFHuge, materials: Materials,
											option: Option) -> None:
	vcf_ref = VCFHuge.read(option.path_ref_VCF)
	iter = chroms_efficients(option)
	first = True
	# とりあえず、後代のVCFも同じ染色体があるとする
	for b, prog_chr, vcf_chr, gmap in zip(iter, vcf.divide_into_chromosomes(),
											vcf_ref.divide_into_chromosomes(),
											materials.chr_maps):
		if not b:
			continue
		vcf_imputed = impute_prog_vcf_chr(vcf_chr, prog_chr, gmap, option)
		with open(option.path_out, 'w' if first else 'a') as out:
			vcf_imputed.write(out, with_header=first)
		first = False

def impute_vcf(option: Option) -> None:
	option.print_info()
	materials = Materials.create(option.path_map, option.path_ped)
	materials.display_map_info()
	vcf = VCFHuge.read(option.path_VCF)
	
	if option.exists_ref():
		impute_progenies(vcf, materials, option)
	else:
		impute_all(vcf, materials, option)


#################### main ####################

option = Option.create(sys.argv)
if option is None:
	Option.usage()
	exit(1)

try:
	impute_vcf(option)
	exit(0)
except ExceptionWithCode as e:
	print(str(e), file=sys.stderr)
	exit(e.get_error_code().value)
