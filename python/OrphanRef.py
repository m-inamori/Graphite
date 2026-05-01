from __future__ import annotations

# coding: utf-8
# OrphanRef.py
# 孤立したサンプルをimputeする

from typing import Optional

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
import Orphan
import RefCommon
from KnownFamily import *
from OptionSmall import OptionSmall


#################### OrphanRef ####################

def impute(samples: list[str], orig_vcf: VCFSmall, ref_haps: list[list[int]],
				phased_vcf: VCFGeno, op: OptionSmall) -> Optional[VCFGenoBase]:
	if not samples:
		return None
	
	vcf = VCFGeno.extract_samples(samples, orig_vcf);
	records = RefCommon.expand_records(vcf, phased_vcf)
	vcf1 = Orphan.impute_samples(samples, records, ref_haps, orig_vcf, op)
	if vcf1 is not None:
		print("%d orphan samples have been imputed." % len(samples))
	return vcf1

__all__ = ['impute']
