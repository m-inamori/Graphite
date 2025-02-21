from __future__ import annotations

# coding: utf-8
# OneKnownFamily.py
# 片親がphasingされていないがknownでもう片親はunknownな家系を補完する

from VCFFamily import *
from Map import *
from KnownFamily import *
from VCFOneParentKnown import VCFOneParentKnown
from VCFOneParentImputed import VCFOneParentImputed

# HMMにrefを使っても計算量が十分に小さいか

def is_small(ref_haps: list[list[int]], L: int) -> bool:
	M = len(ref_haps[0])				# マーカー数
	NH = len(ref_haps)
	R = NH**2 * (2*NH - 1)
	return R * M < 10**8 and R < 10**5 and L * R * M < 10**9

def impute(orig_vcf: VCFSmall, imputed_vcf: VCFSmall, ref_haps: list[list[int]],
						families: list[KnownFamily], gmap: Map) -> VCFSmallBase:
	vcfs: list[VCFSmallBase] = []
	for family in families:
		vcf = VCFFamily.create_by_two_vcfs(imputed_vcf,
											orig_vcf, family.samples())
		# KnownFamilyを使うように直したい
		is_mat_known = family.mat != '0'
		if is_small(ref_haps, len(families)):
			vcf1 = VCFOneParentKnown(vcf.header, vcf.records,
										ref_haps, is_mat_known, gmap)
			vcf1.impute_known_parent()
			for prog in vcf1.samples[2:]:
				samples = [family.mat, family.pat, prog]
				vcf2 = VCFFamily.create_by_two_vcfs(imputed_vcf,
													orig_vcf, samples)
				vcf3 = VCFOneParentImputed(vcf2.header, vcf2.records,
												ref_haps, is_mat_known, gmap)
				vcf3.impute()
				vcfs.append(vcf3)
	
	new_vcf = VCFSmall.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
