from __future__ import annotations

# coding: utf-8
# BothKnownFamilyRef.py
# 両親が分っているがphasingされていない家系を補完する

from functools import reduce
from collections import defaultdict, Counter
from typing import List, Tuple, Optional, IO, Dict, Iterator, Sequence

from VCF import VCFSmall
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import *
import BothKnownFamily
from Map import *
import RefCommon
from OptionSmall import OptionSmall


#################### BothKnownFamilyRef ####################

def impute(orig_vcf: VCFSmall, phased_vcf: VCFGeno, ref_haps: list[list[int]],
			families: list[Family], op: OptionSmall) -> Optional[VCFGenoBase]:
	vcfs: list[VCFGenoBase] = []
	for family in families:
		vcf = VCFGeno.extract_samples(family.samples(), orig_vcf);
		records = RefCommon.merge_family_records(phased_vcf,
													vcf, family.samples())
		family_vcf = BothKnownFamily.impute_family(family, records,
													len(families),
													ref_haps, orig_vcf, op)
		if family_vcf is not None:
			vcfs.append(family_vcf)
	
	if not vcfs:
		return None
	
	print("%d families whose parents are known have been imputed." %
															len(families))
	new_vcf = VCFGeno.join(vcfs, orig_vcf.samples)
	return new_vcf

__all__ = ['impute']
