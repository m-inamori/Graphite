# coding: utf-8
# Matrix
from itertools import count

class Matrix:
	def __init__(self, a):
		self.a	= a
	
	def copy(self):
		return Matrix([ v[:] for v in self.a ])
	
	def __add__(self, B):
		C = [ [ e1 + e2 for e1, e2 in zip(v1, v2) ]
						for v1, v2 in zip(self.a, B.a) ]
		return Matrix(C)
	
	def __sub__(self, B):
		C = [ [ e1 - e2 for e1, e2 in zip(v1, v2) ]
						for v1, v2 in zip(self.a, B.a) ]
		return Matrix(C)
	
	def __mul__(self, B):
		if isinstance(B, Matrix):
			C = [ [ sum(self.a[i][k] * B.a[k][j] for k in range(self.W()))
												 for j in range(B.W()) ]
												 for i in range(self.H()) ]
			return Matrix(C)
		else:	# float or int
			return Matrix([ [ e * B for e in v ] for v in self.a ])
	
	def __truediv__(self, d):
		return Matrix([ [ e / d for e in v ] for v in self.a ])
	
	# e: zero or positive integer
	def __pow__(self, e):
		if e == 0:
			return Matrix([ [ 1 if i == j else 0 for j in range(len(self.a)) ]
												 for i in range(len(self.a)) ])
		elif e == 1:
			return self
		elif e % 2 == 1:
			return self * (self ** (e - 1))
		else:
			M = self ** (e // 2)
			return M * M
	
	def __iadd__(self, B):
		for i in range(self.H()):
			for j in range(self.W()):
				self.a[i][j] += B.a[i][j]
		return self
	
	def __isub__(self, B):
		for i in range(self.H()):
			for j in range(self.W()):
				self.a[i][j] -= B.a[i][j]
		return self
	
	def __imul__(self, B):
		a = [ v[:] for v in self.a ]
		for i in range(self.H()):
			for j in range(self.W()):
				self.a[i][j] = sum(a[i][k] * B.a[k][j] for k in range(self.W()))
		return self
	
	def __str__(self):
		v = [ ' '.join(map(str, v)) for v in self.a ]
		return '\n'.join(v)
	
	def H(self):
		return len(self.a)
	
	def W(self):
		return len(self.a[0])
	
	def is_square(self):
		return self.H() == self.W()
	
	def log(self, limit_time=100):
		assert(self.is_square())
		# log M = log(E - A)
		L = len(self.a)
		A = Matrix.identity(L) - self
		B = Matrix.zero(L, L)	# result
		eps = 1e-16
		P = A.copy()	# A^e
		for e in count(1):
			if e == limit_time:
				break
			Q = P / e
			B -= Q
			if Q.max_abs() < eps:
				break
			P *= A
		
		return B
	
	def exp(self):
		assert(self.is_square())
		L = len(self.a)
		B = Matrix.zero(L, L)	# result
		P = Matrix.identity(L)	# A^e / e!
		eps = 1e-16
		for e in count():
			B += P
			if P.max_abs() < eps:
				break
			P *= self
			P /= e + 1
		return B
	
	def max_abs(self):
		return max(abs(e) for v in self.a for e in v)
	
	@staticmethod
	def identity(L):
		C = [ [ 1.0 if j == i else 0.0 for j in range(L) ] for i in range(L) ]
		return Matrix(C)
	
	@staticmethod
	def zero(H, W):
		return Matrix([ [0.0] * W for _ in range(H) ])
	
	@staticmethod
	def sum(Ms):
		assert(Ms)
		S = Ms[0].copy()
		for i in range(1, len(Ms)):
			S += Ms[i]
		return S
	
	@staticmethod
	def average(Ms):
		return Matrix.sum(Ms) / len(Ms)

