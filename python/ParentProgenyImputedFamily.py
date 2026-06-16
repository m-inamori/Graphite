from __future__ import annotations

# coding: utf-8
# ParentProgenyImputedFamily.py
# 片親と後代がphasingされて片親は分っているがphasingされていない家系を補完する

from typing import Optional, Sequence

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import VCFFamily
from VCFImputable import *
from VCFOneParentImputed import VCFOneParentImputed
from KnownFamily import KnownFamily
from Map import *
from OptionSmall import OptionSmall


#################### ParentProgenyImputedFamily ####################

def impute(orig_vcf: VCFSmall, imputed_vcf: VCFGeno,
									ref_haps: list[list[int]],
									families: list[KnownFamily],
									non_imputed_parents: list[str],
									op: OptionSmall) -> Optional[VCFGenoBase]:
	vcfs: list[VCFImputable] = []
	for family in families:
		vcf_family = VCFFamily.create_by_two_vcfs(imputed_vcf,
													orig_vcf, family.samples())
		should_impute_mat = vcf_family.mat() in non_imputed_parents
		vcf1 = VCFOneParentImputed(family.samples(), vcf_family.records,
													ref_haps, should_impute_mat,
													op.map, orig_vcf)
		vcfs.append(vcf1)
	
	if not vcfs:
		return None
	
	for vcf in vcfs:
		vcf.impute()
	
	print("%d families whose one parent is imputed and"
					" one progeny is imputed have been imputed." % len(vcfs))
	vcfs_base: list[VCFGenoBase] = [ vcf for vcf in vcfs ]
	new_vcf = VCFGeno.join(vcfs_base, orig_vcf.samples)
	return new_vcf

__all__ = ['impute', 'impute_family']
