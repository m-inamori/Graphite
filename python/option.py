from __future__ import annotations

# coding: utf-8
# option.py

from typing import Optional
import sys

#################### Option ####################

class Option:
	def __init__(self, VCF: str, ped: str, m: str,
						families: list[int], chroms: list[int],
						num: int, lp: int, ol: bool, out: str):
		self.path_VCF: str				= VCF
		self.path_ped: str				= ped
		self.path_map: str				= m
		self.families: list[int]		= families
		self.chroms: list[int]			= chroms
		self.num_threads: int			= num
		self.path_out: str				= out
		self.ratio: float				= 0.01
		self.lower_progs: int			= lp
		self.only_large_families: bool	= ol
	
	@staticmethod
	def create(argv: list[str]) -> Optional[Option]:
		def exists(s: str) -> bool:
			return any(argv[i] == s for i in range(4, len(argv) - 1))
		
		def flag_value(s: str) -> str:
			for i in range(3, len(argv) - 2):
				if argv[i] == s:
					return argv[i+1]
			else:
				return ''
		
		# str -> [index]
		def parse_indices(c: str) -> list[int]:
			indices = []
			v = c.split(',')
			for s in v:
				w = s.split(':')
				if len(w) == 1:
					indices.append(int(w[0]))
				elif len(w) == 2:
					indices.extend(list(range(int(w[0]), int(w[1]))))
				else:
					raise ValueError
			return indices
		
		def get_families() -> list[int]:
			s = flag_value('-f')
			if s == '':
				return []	# all
			else:
				return parse_indices(s)
		
		def get_chroms() -> list[int]:
			s = flag_value('-c')
			if s == '':
				return []	# all
			else:
				return parse_indices(s)
		
		def get_num_threads() -> int:
			s = flag_value('-t')
			if s == '':
				return 1
			else:
				return int(s)
		
		def get_lower_progenies() -> int:
			s = flag_value('-p')
			if s == '':
				return 10
			else:
				return int(s)
		
		if not (4 <= len(argv) <= 13):
			return None
		
		try:
			map_path = flag_value('-m')
			families = get_families()
			chroms = get_chroms()
			num_threads = get_num_threads()
			lower_progs = get_lower_progenies()
			ol = exists('-l')
		except ValueError:
			return None
		
		return Option(argv[1], argv[2], map_path, families,
								chroms, num_threads, lower_progs, ol, argv[-1])
	
	@staticmethod
	def usage():
		u = ('python graphite.py VCF ped [-m map] [-t num_threads] ' +
				'[-f family indices] [-c chrom indices] [-l] out.')
		print('usage : ' + u, file=sys.stderr)
		print('indices: (index|first:last)[,(index|first:last)[,..]]',
															file=sys.stderr)


#################### OptionImpute ####################

class OptionImpute:
	def __init__(self, max_dist: int, min_positions: int,
									min_graph: int, min_c: float):
		self.max_dist		= max_dist
		self.min_positions	= min_positions
		self.min_graph		= min_graph
		self.min_crossover: float	= 1.0
