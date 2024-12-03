# coding: utf-8
# error_codes.py

from __future__ import annotations

from Map import Map
from pedigree import PedigreeTable


#################### Materials ####################

class Materials:
	def __init__(self, path_m: str, m: Map, path_p: str, p: PedigreeTable):
		self.path_map: str = path_m
		self.geno_map: Map = m
		self.chr_maps: list[Map] = self.geno_map.create_chr_maps()
		
		self.path_ped: str = path_p
		self.ped: PedigreeTable = p
	
	def iter_chr_maps(self) -> Iterator[Map]:
		if self.chr_maps:
			return iter(self.chr_maps)
		else:
			while True:
				yield Map.default_map()
	
	def total_cM(self) -> float:
		return sum(m.total_cM() for m in self.chr_maps)
	
	def display_map_info(self):
		s = "Genetic Map : "
		if self.geno_map.is_empty():
			s += "default map(1Mbp=1cM)."
		else:
			s += self.path_map + '\n'
			s += "%d chrmosomes %f cM." % (len(self.chr_maps), self.total_cM())
		print(s)
	
	@staticmethod
	def create(path_map: str, path_ped: str) -> Materials:
		geno_map = Map.read(path_map)
		ped = PedigreeTable.read(path_ped)
		return Materials(path_map, geno_map, path_ped, ped)
