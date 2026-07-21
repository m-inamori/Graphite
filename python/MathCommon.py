# coding: utf-8
# MathCommon.py

from __future__ import annotations
from typing import Callable


#################### MathCommon ####################

def Simpson(f: Callable[[float], float], a: float, b: float, n: int) -> float:
	s = f(a) + f(b)
	s += sum(f((k*a+(n*2-k)*b)/(n*2)) for k in range(1, n*2, 2)) * 4
	s += sum(f((k*a+(n*2-k)*b)/(n*2)) for k in range(2, n*2, 2)) * 2
	return s * (b-a) / (n*2) / 3
