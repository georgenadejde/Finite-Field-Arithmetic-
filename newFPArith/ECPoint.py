from Util import *
from ComplexPoint import *
from FP_arith import *


class ECPoint:
	def __init__(self, x, p, y = 0, z = 1):
		
		'''Creation of an elliptic curve point
		Args:
			x (ComplexPoint): x coordinate
			p (int): modulo
			y (ComplexPoint): y coordinate (default: `0`)
			z (ComplexPoint): z coordinate (default: `1`)
		'''

		# y is not needed for montgomery curves
		if y == 0:
			self.Y = ComplexPoint(0,0)
		else:
			self.Y = y

		self.p = p
		self.X = x

		if z == 1:
			self.Z = ComplexPoint(1,0)
		else:
			self.Z = z		
		
	def simplify_point(self):
		# x = X // Z
		self.X = fp2_div(self.X, self.Z, self.p)
		self.Z = ComplexPoint(1, 0)

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
		return self.Z == ComplexPoint(0,0)

	def is_T(self):
		return self.X == ComplexPoint(0,0) and not (self.Z == ComplexPoint(0,0))

	def __str__(self):
		
		return "(X:Z) = ({}:{})".format(self.X,self.Z)
