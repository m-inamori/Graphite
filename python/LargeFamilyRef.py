# coding: utf-8
# LargeFamilyRef.py

from __future__ import annotations
from itertools import *
import sys
from typing import Optional

from VCF import *
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import VCFFamilyRecord, VCFFamilyBase
from VCFImputable import *
from VCFLargeBothRef import *
from VCFLargeOneRef import *
from VCFLargeNoRef import *
import LargeFamily
import RefCommon
from Map import *
from KnownFamily import KnownFamily
from option import *
from Genotype import Genotype
from common import *


#################### process ####################

# 親がphasedかどうかで異なるクラスのVCFにする
# リファレンスにしかないpositionは./.で埋めて
# リファレンスにないpositionは捨てる
def create_family_vcf(vcf_family: VCFFamilyBase,
							ref_vcf: VCFGeno, gmap: Map) -> VCFImputable:
	samples = vcf_family.samples
	records = RefCommon.merge_family_records(ref_vcf, vcf_family, samples)
	is_mat_ref = vcf_family.mat() in ref_vcf.samples
	is_pat_ref = vcf_family.pat() in ref_vcf.samples
	if is_mat_ref and is_pat_ref:
		return VCFLargeBothRef(samples, records, gmap, ref_vcf.vcf)
	elif (not is_mat_ref) and (not is_pat_ref):
		return VCFLargeNoRef(vcf_family.samples, records, gmap, ref_vcf)
	else:
		return VCFLargeOneRef(samples, records, gmap, is_mat_ref, ref_vcf)

def impute(families: list[KnownFamily], orig_vcf: VCFSmall, ref_vcf: VCFGeno,
							gmap: Map, option: Option) -> Optional[VCFGeno]:
	vcfs_family = LargeFamily.impute_all_families(orig_vcf, families,
															gmap, option)
	vcfs: list[VCFImputable] = []
	for vcf_family in vcfs_family:
		vcfs.append(create_family_vcf(vcf_family, ref_vcf, gmap))
	
	for vcf in vcfs:
		vcf.impute()
	
	imputed_vcfs: list[VCFGenoBase] = [ vcf for vcf in vcfs ]
	if imputed_vcfs:
		return VCFGeno.join(imputed_vcfs, orig_vcf.samples)
	else:
		return None

__all__ = ['impute']
