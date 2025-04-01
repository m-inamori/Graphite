# coding: utf-8
# log.py

from functools import reduce
from math import log, exp
from typing import Iterator


#################### Log ####################

LOGZERO = -10000.0

def log_add(log_a: float, log_b: float) -> float:
	# e.g. log_a = log2 log_b = log3
	# log(2 + 3) = log3 + log(2/3 + 1) = log3 + log(exp(log2-log3) + 1)
	if log_a < log_b:
		return log_b + log(exp(log_a - log_b) + 1.0)
	else:
		return log_a + log(exp(log_b - log_a) + 1.0)

def log_sub(log_a: float, log_b: float) -> float:
	# e.g. log_a = log3 log_b = log2
	# log(3 - 2) = log3 + log(1 - 2/3) = log3 + log(1 - exp(log2-log3))
	assert(log_a >= log_b)
	if log_a == log_b or exp(log_b - log_a) == 1.0:
		return LOGZERO
	
	return log_a + log(1.0 - exp(log_b - log_a))

def log_sum(iterable: Iterator[float]) -> float:
	return reduce(log_add, iterable)

def log_diff(log_a: float, log_b: float) -> float:
	if log_a > log_b:
		return log_sub(log_a, log_b)
	else:
		return log_sub(log_b, log_a)

def modified_log(x: float) -> float:
	return LOGZERO if x == 0.0 else log(x)


