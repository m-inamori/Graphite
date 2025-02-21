# coding: utf-8
# pedigree.py

from __future__ import annotations
from typing import Iterator, Optional
import sys

import error_codes
from exception_with_code import *
from common import read_csv


#################### FormatException ####################

class FormatException(ExceptionWithCode):
	def __init__(self, lines: list[str]):
		super().__init__(FormatException.create_message(lines))
	
	def get_error_code(self) -> error_codes.Type:
		return error_codes.Type.PEDIGREE_INVALID_FORMAT
	
	@staticmethod
	def create_message(lines: list[str]) -> str:
		if len(lines) == 1:
			s = "error : the following line doesn't have four columns :"
		else:
			s = "error : the following lines don't have four columns :"
		
		for line in lines:
			s += '\n' + line
		return s


#################### ParentsException ####################

class ParentsException(ExceptionWithCode):
	def __init__(self, lines: list[str]):
		super().__init__(ParentsException.create_message(lines))
	
	def get_error_code(self) -> error_codes.Type:
		return error_codes.Type.PARENT_NOT_DEFINED
	
	@staticmethod
	def create_message(parents: list[str]) -> str:
		if len(parents) == 1:
			s = "error : the following parent isn't defined :"
		else:
			s = "error : the following parents aren't defined :"
		
		for parent in parents:
			s += '\n' + parent
		return s


#################### SamplesException ####################

class SamplesException(ExceptionWithCode):
	def __init__(self, lines: list[str]):
		super().__init__(SamplesException.create_message(lines))
	
	def get_error_code(self) -> error_codes.Type:
		return error_codes.Type.SAMPLES_NOT_IN_PEDIGREE
	
	@staticmethod
	def create_message(samples: list[str]) -> str:
		if len(samples) == 1:
			s = "error : the following sample isn't in pedigree :"
		else:
			s = "error : the following samples aren't in pedigree :"
		
		for sample in samples:
			s += '\n' + sample
		return s


#################### Family ####################

class Family:
	def __init__(self, mat: str, pat: str, progs: list[str]):
		self.mat: str				= mat
		self.pat: str				= pat
		self.progenies: list[str]	= progs
	
	def parents(self) -> tuple[str, str]:
		return (self.mat, self.pat)
	
	def num_progenies(self) -> int:
		return len(self.progenies)
	
	def samples(self) -> list[str]:
		return [self.mat, self.pat] + self.progenies


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
						if p.mat != '0' and p.pat != '0' and
							p.mat in samples and p.pat in samples ]
		return PedigreeTable(progs)
	
	def generate_parents(self) -> Iterator[tuple[str, str]]:
		return (p.parents() for p in self.table)
	
	def get_children(self, parents: tuple[str, str]) -> list[str]:
		return [ p.name for p in self.table if p.parents() == parents ]
	
	def extract_families(self) -> list[Family]:
		parents = sorted(set((prog.mat, prog.pat) for prog in self.table))
		return [ Family(*p, self.get_children(p)) for p in parents ]
	
	def limit_samples(self, samples: list[str]) -> PedigreeTable:
		ped_samples = set(p.name for p in self.table)
		missing_samples = [ s for s in samples if s not in ped_samples ]
		if missing_samples:
			raise SamplesException(missing_samples)
		
		set_samples = set(samples)
		
		progs = [ p for p in self.table if p.name in set_samples ]
		return PedigreeTable(progs)
	
	def check_parents(self) -> list[str]:
		missing_parents = set()
		set_progs = set(prog.name for prog in self.table)
		for prog in self.table:
			for p in prog.parents():
				if p != '0' and p not in set_progs:
					missing_parents.add(p)
		return list(missing_parents)
	
	@staticmethod
	def read_lines(path: str) -> list[list[str]]:
		table: list[list[str]] = []
		not_four_columns_lines: list[str] = []
		try:
			with open(path, 'r') as f:
				for line in f:
					v = line.split()
					if len(v) != 4:
						not_four_columns_lines.append(line.rstrip())
					else:
						table.append(v)
			if not_four_columns_lines:
				raise FormatException(not_four_columns_lines)
		except FileNotFoundError:
			raise FileNotFoundException(path)
		except IOError:
			raise FileNotFoundException(path)
		
		return table
	
	@staticmethod
	def remove_duplidated_progenies(progs: list[Progeny]) -> list[Progeny]:
		new_progs: list[Progeny] = []
		s = set()
		for prog in progs:
			if prog.name not in s:
				new_progs.append(prog)
				s.add(prog.name)
			else:
				print('Progeny %s is duplicated.' % prog.name)
		return new_progs
	
	@staticmethod
	def create(progs: list[Progeny]) -> PedigreeTable:
		ped = PedigreeTable(progs)
		missing_parents = ped.check_parents()
		if missing_parents:
			raise ParentsException(missing_parents)
		return ped
	
	@staticmethod
	def read(path: str) -> PedigreeTable:
		table = PedigreeTable.read_lines(path)
		progs: list[Progeny] = [ Progeny(v) for v in table ]
		new_progs = PedigreeTable.remove_duplidated_progenies(progs)
		return PedigreeTable.create(new_progs)
