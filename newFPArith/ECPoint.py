from Util import *
from Fp2Point import *
from FP_arith import *


class ECPoint:
	def __init__(self, x, p, y = 0, z = 1):
		
		'''Creation of an elliptic curve point
		Args:
			x (Fp2Point): x coordinate
			p (int): modulo
			y (Fp2Point): y coordinate (default: `0`)
			z (Fp2Point): z coordinate (default: `1`)
		'''

		# y is not needed for montgomery curves
		if y == 0:
			self.Y = Fp2Point(0,0)
		else:
			self.Y = y

		self.p = p
		self.X = x

		if z == 1:
			self.Z = Fp2Point(1,0)
		else:
			self.Z = z		
		
	def normalize(self):
		# x = X // Z
		self.X = fp2_div(self.X, self.Z, self.p)
		self.Z = Fp2Point(1, 0)

		return self

	def __eq__(self, other):
		if isinstance(other, ECPoint):
 			return self.X == other.X and self.Z == other.Z
		return False	

	def X(self):
		return self.X

	def Z(self):
		return self.Z

	def p(self):
		return self.p

	def is_POIF(self):
		return self.Z == Fp2Point(0,0)

	def is_T(self):
		return self.X == Fp2Point(0,0) and not (self.Z == Fp2Point(0,0))

	def __str__(self):
		
		return "(X:Z) = ({}:{})".format(self.X,self.Z)
