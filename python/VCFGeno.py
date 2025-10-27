from __future__ import annotations

# coding: utf-8
# VCFGeno.py
# 基本的にGenotypeのみ持つ

from itertools import count
from abc import ABC, abstractmethod
from typing import TextIO

from GenoRecord import GenoRecord
from VCF import VCFSmall, VCFRecord
from Genotype import Genotype


#################### VCFGenoBase ####################

class VCFGenoBase(ABC):
	def __init__(self, samples: list[str], vcf: VCFSmall) -> None:
		self.samples: list[str] = samples
		self.vcf: VCFSmall = vcf
	
	@abstractmethod
	def __len__(self) -> int:
		pass
	
	@abstractmethod
	def get_record(self, i: int) -> GenoRecord:
		pass
	
	def get_ref_vcf(self) -> VCFSmall:
		return self.vcf
	
	def get_samples(self) -> list[str]:
		return self.samples
	
	def num_samples(self) -> int:
		return len(self.samples)
	
	def extract_columns(self, samples: list[str]) -> list[int]:
		dic = dict(zip(self.samples, count()))
		return [ dic.get(sample, -1) for sample in samples ]
	
	def clip_raw_haplotype(self, sample_index: int, side: int) -> list[int]:
		hap: list[int] = []
		for j in range(len(self)):
			record = self.get_record(j)
			h = record.get_allele(sample_index, side)
			hap.append(h)
		return hap
	
	def write(self, out: TextIO, with_header: bool=True) -> None:
		if with_header:
			self.vcf.write_header(out)
		for i in range(len(self)):
			self.get_record(i).write(self.vcf.records[i], out)
	
	def extract_by_samples(self, samples: list[str]) -> VCFGeno:
		cs = self.extract_columns(samples)
		new_records: list[GenoRecord] = []
		for i in range(len(self)):
			record = self.get_record(i)
			pos = record.pos
			geno = [ record.geno[c] for c in cs ]
			new_record = GenoRecord(pos, geno)
			new_records.append(new_record)
		return VCFGeno(samples, new_records, self.vcf)


#################### VCFGeno ####################

class VCFGeno(VCFGenoBase):
	def __init__(self, samples: list[str],
						records: list[GenoRecord], vcf: VCFSmall) -> None:
		super().__init__(samples, vcf)
		self.records: list[GenoRecord] = records
	
	def __len__(self) -> int:
		return len(self.records)
	
	def get_record(self, i: int) -> GenoRecord:
		return self.records[i]
	
	@staticmethod
	def extract_samples(samples: list[str], vcf: VCFSmall) -> VCFGeno:
		cs = vcf.extract_columns(samples)
		new_records: list[GenoRecord] = []
		for i in range(len(vcf)):
			record = vcf.get_record(i)
			pos = int(record.v[1])
			new_v = [ Genotype.all_gt_to_int(record.v[c]) for c in cs ]
			new_record = GenoRecord(pos, new_v)
			new_records.append(new_record)
		return VCFGeno(samples, new_records, vcf)
	
	# samplesの順番で結合
	@staticmethod
	def join(vcfs: list[VCFGenoBase], samples: list[str]) -> VCFGeno:
		dic = { }
		for vcf in vcfs:
			# sに'0'があっても問題ない
			for c, s in enumerate(vcf.samples):
				dic[s] = (vcf, c)
		cols = []
		for s in samples:
			if s in dic:
				cols.append((s,) + dic[s])
		
		new_samples = [ s for s, vcf, c in cols ]
		new_records = []
		for i in range(len(vcfs[0])):
			geno = []
			for s, vcf, c in cols:
				geno.append(vcf.get_record(i).geno[c])
			new_records.append(GenoRecord(vcfs[0].get_record(i).pos, geno))
		return VCFGeno(new_samples, new_records, vcfs[0].vcf)
	
	# vcf1にあるsampleはvcf1から、
	# そうでないsampleはvcf2からGenotypeを取って新たなVCFを作る
	@staticmethod
	def create_by_two_vcfs(vcf1: VCFGenoBase, vcf2: VCFSmall,
										samples: list[str]) -> VCFGeno:
		# samplesは[mat, pat, prog1, ...]の前提
		columns1 = vcf1.extract_columns(samples)
		columns2 = vcf2.extract_columns(samples)
		new_records: list[GenoRecord] = []
		for i in range(len(vcf1)):
			record1: GenoRecord = vcf1.get_record(i)
			record2: VCFRecord = vcf2.get_record(i)
			geno: list[int] = []
			for j in range(len(samples)):
				if columns1[j] != -1:
					geno.append(record1.geno[columns1[j]])
				elif columns2[j] != -1:
					geno.append(Genotype.all_gt_to_int(record2.v[columns2[j]]))
				else:	# どちらのVCFにもそのサンプルは無い
					geno.append(Genotype.NA)
			new_record = GenoRecord(record1.pos, geno)
			new_records.append(new_record)
		return VCFGeno(samples, new_records, vcf2)
