from __future__ import annotations

# coding: utf-8
# OptionSmall.py

from Map import *


#################### OptionSmall ####################

class OptionSmall:
	def __init__(self, geno_map: Map, num_t: int, pratio: float,
												ii: bool, ou: bool) -> None:
		self.map: Map = geno_map
		self.num_threads: int = num_t
		self.precision_ratio: float = pratio
		self.imputes_isolated_samples: bool = ii
		self.outputs_unimputed_samples: bool = ou
