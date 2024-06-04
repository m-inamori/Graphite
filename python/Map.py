from __future__ import annotations

# coding: utf-8
# Map.py

from itertools import *
from typing import Iterator, List, Any, IO

import error_codes
from exception_with_code import ExceptionWithCode, FileNotFoundException
from common import *


#################### Map ####################

class Map:
	class Record:
		def __init__(self, chr: str, cM: float, Mbp: float):
			self.chr: str	= chr
			self.cM: float	= cM
			self.Mbp: float	= Mbp
		
		def __str__(self) -> str:
			return "%s,%f,%f" % (self.chr, self.cM, self.Mbp)
		
		@staticmethod
		def create(v: list[str]) -> Map.Record:
			try:
				record = Map.Record(v[0], float(v[1]), float(v[2]))
				return record
			except:
				lines = [ '.'.join(v) ]
				raise MapFormatException(lines)
	
	def __init__(self, records: List['Record']):
		self.records: List['Map.Record'] = records
	
	def is_empty(self) -> bool:
		return not self.records
	
	def total_cM(self) -> float:
		return self.records[-1].cM
	
	def bp_to_cM(self, bp: float) -> float:
		Mbp: float = bp * 1e-6
		if self.is_empty():
			return Mbp
		
		index: int = self.binary_search(0, len(self.records)-1, Mbp)
		r1: Map.Record = self.records[index]
		r2: Map.Record = self.records[index+1]
		if r1.Mbp == r2.Mbp:
			return r1.cM
		return r1.cM + (r2.cM-r1.cM) / (r2.Mbp-r1.Mbp) * (Mbp-r1.Mbp)
	
	def binary_search(self, first: int, last: int, Mbp: float) -> int:
		if first == last - 1:
			return first
		
		mid: int = (first + last) // 2
		if Mbp < self.records[mid].Mbp:
			return self.binary_search(first, mid, Mbp)
		else:
			return self.binary_search(mid, last, Mbp)
	
	def create_chr_maps(self) -> list[Map]:
		if self.is_empty():
			# Make one map and use it for every chromosome.
			# That makes it easier to delete.
			return []
		else:
			return list(self.divide_into_chromosomes())
	
	def check_in_order(self):
		past_chrs: list[str] = []
		for r1, r2 in zip(self.records, self.records[1:]):
			if r1.chr != r2.chr:
				past_chrs.append(r1.chr)
				if r2.chr in past_chrs:
					raise TwiceChrException(r2.chr)
			else:
				if r1.cM > r2.cM:
					raise OutOfOrderCMException(r1, r2)
				elif r1.Mbp > r2.Mbp:
					raise OutOfOrderMbpException(r1, r2)
	
	def divide_into_chromosomes(self) -> Iterator[Map]:
		for chr, v in groupby(self.records, key=lambda r: r.chr):
			yield Map(list(v))
	
	@staticmethod
	def read_lines(path: str) -> list[list[str]]:
		try:
			table: list[list[str]] = []
			with open(path) as f:
				not_three_columns_lines: list[str] = []
				for line in f:
					v = line.split(',')
					if len(v) != 3 and len(not_three_columns_lines) < 5:
						not_three_columns_lines.append(line)
					table.append(v)
			if not_three_columns_lines:
				raise MapFormatException(not_three_columns_lines)
		except FileNotFoundError:
			raise FileNotFoundException(path)
		except IOError:
			raise FileNotFoundException(path)
		
		return table
	
	@staticmethod
	def read(path: str) -> Map:
		if path:
			table = Map.read_lines(path)
			records = [ Map.Record.create(v) for v in table ]
		else:
			records = []
		geno_map = Map(records)
		geno_map.check_in_order()
		return geno_map
	
	@staticmethod
	def default_map() -> Map:
		r1 = Map.Record('1', 0.0, 0.0)
		r2 = Map.Record('1', 1.0, 1.0)
		return Map([r1, r2])


#################### VCFMeasurable ####################

class VCFMeasurable(object):
	def __init__(self, map_: Map):
		self.map: Map = map_
	
	def cM(self, pos: int) -> float:
		return self.map.bp_to_cM(pos)


#################### MapFormatException ####################

class MapFormatException(ExceptionWithCode):
	def __init__(self, lines: list[str]):
		super().__init__(MapFormatException.create_message(lines))
	
	def get_error_code(self) -> error_codes.Type:
		return error_codes.Type.MAP_INVALID_FORMAT
	
	@staticmethod
	def create_message(lines: list[str]):
		if len(lines) == 1:
			s = "error : the following line doesn't have three columns :"
		else:
			s = "error : the following lines don't have three columns :"
		
		for line in lines:
			s += '\n' + line
		return s


#################### TwiceChrException ####################

class TwiceChrException(ExceptionWithCode):
	def __init__(self, chr: str):
		super().__init__(TwiceChrException.create_message(chr))
	
	def get_error_code(self) -> error_codes.Type:
		return error_codes.Type.MAP_INVALID_FORMAT
	
	@staticmethod
	def create_message(chr: str):
		s = "error : the following chromosome comes out twice :\n"
		s += chr
		return s


#################### OutOfOrderCMException ####################

class OutOfOrderCMException(ExceptionWithCode):
	def __init__(self, r1: Map.Record, r2: Map.Record):
		super().__init__(OutOfOrderCMException.create_message(r1, r2))
	
	def get_error_code(self) -> error_codes.Type:
		return error_codes.Type.MAP_INVALID_FORMAT
	
	@staticmethod
	def create_message(r1: Map.Record, r2: Map.Record):
		s = "error : the following cMs are out of order :"
		s += str(r1)
		s += str(r2)
		return s


#################### OutOfOrderMbpException ####################

class OutOfOrderMbpException(ExceptionWithCode):
	def __init__(self, r1: Map.Record, r2: Map.Record):
		super().__init__(OutOfOrderMbpException.create_message(r1, r2))
	
	def get_error_code(self) -> error_codes.Type:
		return error_codes.Type.MAP_INVALID_FORMAT
	
	@staticmethod
	def create_message(r1: Map.Record, r2: Map.Record):
		s = "error : the following Mbps are out of order :\n"
		s += str(r1) + '\n' + str(r2)
		return s
