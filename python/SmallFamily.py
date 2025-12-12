from __future__ import annotations

# coding: utf-8
# SmallFamily.py

from itertools import *
from collections import defaultdict, Counter
import random
import csv
import sys
from multiprocessing import Pool

from typing import Dict, List, Tuple, Set, Optional

from VCF import *
from VCFGeno import VCFGenoBase, VCFGeno
from VCFFamily import VCFFamily, VCFFamilyRecord
from pedigree import PedigreeTable, Family
from KnownFamily import KnownFamily
import BothImputedFamily
import BothKnownFamily
import OneKnownFamily
import OneImputedFamily
import ProgenyImputedFamily
import SelfFamily
import SelfNonImputedFamily
import Orphan
import ImputedAndKnownFamily
from VCFIsolated import VCFIsolated
import ReferenceHaplotype
from SampleManager import *
from Map import *
from OptionSmall import OptionSmall
from common import *


#################### SmallFamily ####################

def impute_vcf_by_both_imputed_parents(orig_vcf: VCFSmall,
									imputed_vcf: VCFGenoBase,
									sample_man: SampleManager,
									op_small: OptionSmall) -> Optional[VCFGeno]:
	families = sample_man.extract_both_imputed_families()
	vcf = BothImputedFamily.impute(orig_vcf, imputed_vcf, families, op_small)
	if vcf is None:
		return None
	
	merged_vcf = VCFGeno.join([imputed_vcf, vcf], orig_vcf.samples)
	sample_man.set(vcf.samples)
	return merged_vcf

def impute_vcf_by_imputed_and_known_parent(
							orig_vcf: VCFSmall, imputed_vcf: VCFGeno,
							ref_haps: list[list[int]], op_small: OptionSmall,
							sample_man: SampleManager) -> Optional[VCFGeno]:
	families = sample_man.extract_imputed_and_known_families()
	# families have been already selected
	# in which one parent has been imputed and one parent has not been imputed
	# collect not phased parents
	samples: list[str] = [ parent for family in families
								  for parent in family.parents()
								  if not sample_man.is_imputed(parent) ]
	
	vcf = ImputedAndKnownFamily.impute(orig_vcf, imputed_vcf,
										ref_haps, families, samples, op_small)
	if vcf is None:
		return None
	
	merged_vcf = VCFGeno.join([imputed_vcf, vcf], orig_vcf.samples)
	sample_man.set(vcf.get_samples())
	return merged_vcf

def impute_vcf_by_both_known_parents(orig_vcf: VCFSmall, imputed_vcf: VCFGeno,
							ref_haps: list[list[int]], op_small: OptionSmall,
							sample_man: SampleManager) -> Optional[VCFGeno]:
	families = sample_man.extract_both_known_families()
	vcf = BothKnownFamily.impute(orig_vcf, imputed_vcf,
									ref_haps, families, op_small)
	if vcf is None:
		return None
	
	merged_vcf = VCFGeno.join([imputed_vcf, vcf], orig_vcf.samples)
	sample_man.set(vcf.get_samples())
	return merged_vcf

def impute_vcf_by_imputed_parent(orig_vcf: VCFSmall, imputed_vcf: VCFGeno,
							ref_haps: list[list[int]], op_small: OptionSmall,
							sample_man: SampleManager) -> Optional[VCFGeno]:
	families = sample_man.extract_one_imputed_families()
	vcf = OneImputedFamily.impute(orig_vcf, imputed_vcf,
											ref_haps, families, op_small)
	if vcf is None:
		return None
	
	merged_vcf = VCFGeno.join([imputed_vcf, vcf], orig_vcf.samples)
	sample_man.set(vcf.get_samples())
	return merged_vcf

def impute_vcf_by_known_parent(orig_vcf: VCFSmall, imputed_vcf: VCFGeno,
							ref_haps: list[list[int]], op_small: OptionSmall,
							sample_man: SampleManager) -> Optional[VCFGeno]:
	families = sample_man.extract_one_known_parent_families()
	vcf = OneKnownFamily.impute(orig_vcf, imputed_vcf,
											ref_haps, families, op_small)
	if vcf is None:
		return None
	
	merged_vcf = VCFGeno.join([imputed_vcf, vcf], orig_vcf.samples)
	sample_man.set(vcf.get_samples())
	return merged_vcf

def impute_self_vcf(orig_vcf: VCFSmall, imputed_vcf: VCFGeno,
							ref_haps: list[list[int]], op_small: OptionSmall,
							sample_man: SampleManager) -> Optional[VCFGeno]:
	families = sample_man.extract_small_self_families()
	imputed_samples = [ s for f in families for s in f.samples()
												if sample_man.is_imputed(s) ]
	vcf = SelfFamily.impute(orig_vcf, imputed_vcf, ref_haps,
									families, imputed_samples, op_small)
	if vcf is None:
		return None
	
	merged_vcf = VCFGeno.join([imputed_vcf, vcf], orig_vcf.samples)
	sample_man.set(vcf.get_samples())
	return merged_vcf

def impute_self_non_imputed_vcf(orig_vcf: VCFSmall, imputed_vcf: VCFGeno,
							ref_haps: list[list[int]], op_small: OptionSmall,
							sample_man: SampleManager) -> Optional[VCFGeno]:
	families = sample_man.extract_self_non_imputed_families()
	vcf = SelfNonImputedFamily.impute(orig_vcf, imputed_vcf, ref_haps,
														families, op_small)
	if vcf is None:
		return None
	
	merged_vcf = VCFGeno.join([imputed_vcf, vcf], orig_vcf.samples)
	sample_man.set(vcf.get_samples())
	return merged_vcf

def impute_vcf_by_progenies(orig_vcf: VCFSmall, imputed_vcf: VCFGeno,
							ref_haps: list[list[int]], op_small: OptionSmall,
							sample_man: SampleManager) -> Optional[VCFGeno]:
	families = sample_man.extract_progenies_imputed_families()
	imputed_progenies = [ [ prog for prog in family.progenies
									if sample_man.is_imputed(prog) ]
												for family in families ]
	vcf = ProgenyImputedFamily.impute(orig_vcf, imputed_vcf, families,
										imputed_progenies, ref_haps, op_small)
	if vcf is None:
		return None
	
	merged_vcf = VCFGeno.join([imputed_vcf, vcf], orig_vcf.samples)
	sample_man.set(vcf.samples)
	return merged_vcf

def impute_orphan_samples(orig_vcf: VCFSmall, imputed_vcf: VCFGeno,
							ref_haps: list[list[int]], op_small: OptionSmall,
							sample_man: SampleManager) -> Optional[VCFGeno]:
	samples = sample_man.extract_non_imputed_samples()
	vcf = Orphan.impute(samples, orig_vcf, ref_haps, op_small)
	if vcf is None:
		return None
	
	merged_vcf = VCFGeno.join([imputed_vcf, vcf], orig_vcf.samples)
	sample_man.set(vcf.get_samples())
	return merged_vcf

def impute_non_imputed_samples(orig_vcf: VCFSmall, merged_vcf: VCFGeno,
											sample_man: SampleManager,
											op: OptionSmall) -> VCFGeno:
	samples = sample_man.extract_non_imputed_samples()
	if samples and op.imputes_isolated_samples:
		reference = sample_man.collect_reference()
		vcf = VCFIsolated.create(orig_vcf, merged_vcf, samples, reference, op)
		vcf.impute()
		new_vcf = VCFGeno.join([merged_vcf, vcf], orig_vcf.samples)
		return new_vcf
	elif samples and op.outputs_unimputed_samples:
		vcf_isolated = VCFGeno.extract_samples(samples, orig_vcf)
		new_vcf = VCFGeno.join([merged_vcf, vcf_isolated], orig_vcf.samples)
		return new_vcf
	else:
		return merged_vcf

def impute(orig_vcf: VCFSmall, merged_vcf: VCFGeno,
							op_small: OptionSmall, sample_man: SampleManager,
							imputes_isolated_samples: bool) -> VCFGeno:
	ref_haps = ReferenceHaplotype.extract_haplotypes(merged_vcf, sample_man)
	while True:
		# 両親が補完されているが子どもが少ない家系を補完する
		new_merged_vcf1 = impute_vcf_by_both_imputed_parents(orig_vcf,
														merged_vcf,
														sample_man, op_small)
		if new_merged_vcf1 is not None:
			merged_vcf = new_merged_vcf1
			continue
		
		# 片親が補完されている家系を補完する
		new_merged_vcf2 = impute_vcf_by_imputed_and_known_parent(orig_vcf,
														merged_vcf, ref_haps,
														op_small, sample_man)
		if new_merged_vcf2 is not None:
			merged_vcf = new_merged_vcf2
			continue
		
		# 両親が補完されていない家系を補完する
		new_merged_vcf3 = impute_vcf_by_both_known_parents(orig_vcf,
														merged_vcf, ref_haps,
														op_small, sample_man)
		if new_merged_vcf3 is not None:
			merged_vcf = new_merged_vcf3
			continue
		
		new_merged_vcf4 = impute_vcf_by_imputed_parent(orig_vcf,
														merged_vcf, ref_haps,
														op_small, sample_man)
		if new_merged_vcf4 is not None:
			merged_vcf = new_merged_vcf4
			continue
		
		# Impute families whose progenies have been imputed
		new_merged_vcf5 = impute_vcf_by_progenies(orig_vcf,
													merged_vcf, ref_haps,
													op_small, sample_man)
		if new_merged_vcf5 is not None:
			merged_vcf = new_merged_vcf5
			continue
		
		# 片親がknownだが補完されていなくてもう片親がunknowな家系
		new_merged_vcf6 = impute_vcf_by_known_parent(orig_vcf,
													merged_vcf, ref_haps,
													op_small, sample_man)
		if new_merged_vcf6 is not None:
			merged_vcf = new_merged_vcf6
			continue
		
		# 自殖で一つでもサンプルがimputedな家系
		new_merged_vcf7 = impute_self_vcf(orig_vcf, merged_vcf, ref_haps,
														op_small, sample_man)
		if new_merged_vcf7 is not None:
			merged_vcf = new_merged_vcf7
			continue
		
		# 自殖でimputedなサンプルが一つもない家系
		new_merged_vcf8 = impute_self_non_imputed_vcf(orig_vcf, merged_vcf,
												ref_haps, op_small, sample_man)
		if new_merged_vcf8 is not None:
			merged_vcf = new_merged_vcf8
			continue
		
		if imputes_isolated_samples:
			new_merged_vcf9 = impute_orphan_samples(orig_vcf,
													merged_vcf, ref_haps,
													op_small, sample_man)
			if new_merged_vcf9 is not None:
				merged_vcf = new_merged_vcf9
				continue
		
		break
	
	merged_vcf = impute_non_imputed_samples(orig_vcf, merged_vcf,
													sample_man, op_small)
	return merged_vcf

__all__ = ['impute_small_family_VCFs']
