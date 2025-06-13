from __future__ import annotations

# coding: utf-8
# option.py

from typing import Optional
import sys


#################### Option ####################

class Option:
	def __init__(self, VCF: str, ref: str, ped: str, m: str,
						families: list[int], chroms: list[int],
						num: int, lp: int, pr: float,
						ii: bool, ou: bool, out: str) -> None:
		self.path_VCF: str						= VCF
		self.path_ped: str						= ped
		self.path_ref_VCF: str					= ref
		self.path_map: str						= m
		self.families: list[int]				= families
		self.chroms: list[int]					= chroms
		self.num_threads: int					= num
		self.path_out: str						= out
		self.ratio: float						= 0.01
		self.lower_progs: int					= lp
		self.precision_ratio: float				= pr
		self.imputes_isolated_samples: bool		= ii
		self.outputs_unimputed_samples: bool	= ou
	
	def exists_ref(self) -> bool:
		return self.path_ref_VCF != ''
	
	def print_info(self) -> None:
		# required
		print("input VCF : %s" % self.path_VCF, file=sys.stderr)
		print("pedigree : %s" % self.path_ped, file=sys.stderr)
		print("output VCF : %s" % self.path_out, file=sys.stderr)
		
		# optional
		print("number of threads : %s" % self.num_threads, file=sys.stderr)
		print("number of progenies for large family : %s" % self.lower_progs,
															file=sys.stderr)
		if self.path_ref_VCF:
			print("ref VCF : %s" % self.path_ref_VCF, file=sys.stderr)
		else:
			print("ref VCF is not specified.", file=sys.stderr)
		
		if not self.imputes_isolated_samples:
			print("isolate samples will not be imputed.", file=sys.stderr)
			if self.outputs_unimputed_samples:
				print("but, outputs these samples", file=sys.stderr)
	
	@staticmethod
	def exists(s: str, argv: list[str]) -> bool:
		return s in argv[1:]
	
	@staticmethod
	def flag_value(flag: str, argv: list[str]) -> str:
		for option, value in zip(argv, argv[1:]):
			if option == flag:
				return value
		else:
			return ''
	
	@staticmethod
	def flags_value(flags: list[str], argv: list[str]) -> str:
		for option, value in zip(argv, argv[1:]):
			if any(option == f for f in flags):
				return value
		else:
			return ''
	
	# str -> [index]
	@staticmethod
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
	
	@staticmethod
	def get_families(argv: list[str]) -> list[int]:
		s = Option.flags_value(['-f', '--family'], argv)
		if s == '':
			return []	# all
		else:
			return Option.parse_indices(s)
	
	@staticmethod
	def get_chroms(argv: list[str]) -> list[int]:
		s = Option.flags_value(['-c', '--chrom'], argv)
		if s == '':
			return []	# all
		else:
			return Option.parse_indices(s)
	
	@staticmethod
	def get_num_threads(argv: list[str]) -> int:
		s = Option.flags_value(['-t', '--num-threads'], argv)
		if s == '':
			return 1
		else:
			return int(s)
	
	@staticmethod
	def get_lower_progenies(argv: list[str]) -> int:
		s = Option.flag_value('--lower-progs', argv)
		if s == '':
			return 10
		else:
			return int(s)
	
	@staticmethod
	def get_precision_ratio(argv: list[str]) -> float:
		if Option.exists('--fast', argv):
			return 0.1
		elif Option.exists('--precision', argv):
			return 10.0
		
		s = Option.flag_value('--precision-ratio', argv)
		if s == '':
			return 1.0
		else:
			return float(s)
	
	@staticmethod
	def create(argv: list[str]) -> Optional[Option]:
		try:
			# The top three are required arguments
			path_vcf = Option.flags_value(['-i', '--input'], argv)
			if not path_vcf:
				raise Exception('input VCF not specified.')
			
			path_ped = Option.flags_value(['-p', '--pedigree'], argv)
			if not path_ped:
				raise Exception('pedigree file not specified.')
			
			path_out = Option.flags_value(['-o', '--output'], argv)
			if not path_out:
				raise Exception('output VCF not specified.')
			
			# Optional
			path_map = Option.flags_value(['-m', '--map'], argv)
			ref_vcf = Option.flags_value(['-r', '--ref'], argv)
			families = Option.get_families(argv)
			chroms = Option.get_chroms(argv)
			num_threads = Option.get_num_threads(argv)
			lower_progs = Option.get_lower_progenies(argv)
			prec_ratio = Option.get_precision_ratio(argv)
			impute_isolated = not Option.exists('--not-impute-isolated', argv)
			out_isolated = Option.exists('--out-isolated',  argv)
			if impute_isolated and out_isolated:
				return None
			
			return Option(path_vcf, ref_vcf, path_ped, path_map,
							families, chroms, num_threads, lower_progs,
							prec_ratio, impute_isolated, out_isolated, path_out)
		except ValueError:
			return None
		except Exception as e:
			print('error : ' + str(e), file=sys.stderr)
			return None
	
	@staticmethod
	def usage() -> None:
		messages: list[str] = [
			'Usage:',
			'  python graphite.py [options]',
			'',
			'  -i <path>/--input <path>     Input VCF file',
			'  -p <path>/--pedigree <path>  Pedigree file',
			'  -o <path>/--output <path>    Output file',
			'',
			'Options:',
			'  -m <path>/--map <path>       Input map file',
			'  -r <path>/--ref <path>       Reference VCF file',
			'  -t <int>/--num-threads <int> Number of threads (1)',
			'  --lower-progs <int>          Lower number of progenies',
			'                               considered to a large family',
			'  --not-impute-isolated        Samples that are isolated in the',
			'                               pedigree are not imputed',
			'  --out-isolated               Do not output isolated samples',
			'                               only valid when used in conjunction',
			'                               with --not-impute-isolated',
			'  -c <int>/--chrom <int>       Impute only the chromosome represented',
			'                               by the specified index (0-based).',
			'  --precision-ratio <float>    Control the runtime for small pedigree',
			'                               HMM analysis. Default is 1.0.',
			'                               Larger values increase runtime.',
			'  --fast                       Shortcut for setting --precision-ratio',
			'                               to 0.1. Optimized for faster runtime',
			'                               with reduced precision.',
			'  --precision                  Shortcut for setting --precision-ratio',
			'                               to 10.0. Enhanced precision at the',
			'                               cost of increased runtime.'
		]
		for msg in messages:
			print(msg, file=sys.stderr)


#################### OptionImpute ####################

class OptionImpute:
	def __init__(self, max_dist: int, min_positions: int,
									min_graph: int, min_c: float):
		self.max_dist		= max_dist
		self.min_positions	= min_positions
		self.min_graph		= min_graph
		self.min_crossover: float	= 1.0
