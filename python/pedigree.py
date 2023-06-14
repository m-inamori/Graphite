# coding: utf-8
# pedigree.py

from __future__ import annotations
from typing import Optional, Iterator

from common import read_csv


#################### Family ####################

class Family:
	def __init__(self, mat: str, pat: str, progs: list[str]):
		self.mat: str	= mat
		self.pat: str	= pat
		self.progenies: list[str]	= progs
	
	def parents(self) -> tuple[str, str]:
		return (self.mat, self.pat)
	
	def num_progenies(self) -> int:
		return len(self.progenies)
	
	def samples(self) -> list[str]:
		return [self.mat, self.pat] + self.progenies
	
	def __hash__(self):
		return hash((self.mat, self.pat))


#################### Progeny ####################

class Progeny:
	def __init__(self, v: list[str]):
		self.family: str	= v[0]
		self.name: str		= v[1]
		self.mat: str		= v[2]
		self.pat: str		= v[3]
	
	def parents(self) -> tuple[str, str]:
		return (self.mat, self.pat)


#################### PedigreeTable ####################

class PedigreeTable:
	def __init__(self, progs: list[Progeny]):
		self.table: list[Progeny]	= progs
	
	def __len__(self) -> int:
		return len(self.table)
	
	def filter_with_parents(self, samples: list[str]) -> PedigreeTable:
		progs = [ p for p in self.table
						if p.name in samples and
							p.mat in samples and p.pat in samples ]
		return PedigreeTable(progs)
	
	def generate_parents(self):
		return (p.parents() for p in self.table)
	
	def get_children(self, parents):
		return [ p.name for p in self.table if p.parents() == parents ]
	
	def extract_families(self) -> list[Family]:
		parents = sorted(set((prog.mat, prog.pat) for prog in self.table))
		return [ Family(*p, self.get_children(p)) for p in parents ]
	
	def get_family(self, mat: str, pat: str) -> Optional[Family]:
		progs = [ prog.name for prog in self.table
								if prog.parents() == (mat, pat) ]
		if not progs:
			return None
		
		return Family(mat, pat, progs)
	
	def limit_samples(self, samples: list[str]) -> PedigreeTable:
		# C++にするときは、progをコピーする
		set_samples = set(samples)
		progs = [ prog for prog in self.table if prog.name in set_samples ]
		return PedigreeTable(progs)
	
	@staticmethod
	def read(path: str) -> PedigreeTable:
		progs: list[Progeny] = []
		samples: set[str] = set()
		for v in read_csv(path, ' '):
			if v[1] not in samples:
				progs.append(Progeny(v))
				samples.add(v[1])
		
		return PedigreeTable(progs)
