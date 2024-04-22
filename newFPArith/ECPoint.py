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
			self.Y = ComplexPoint(0,0,p)
		else:
			self.Y = y

		self.p = p

		# poif in projective coordinates
		if not x and not z:
			self.X = ComplexPoint(1, 0, p)
			self.Z = ComplexPoint(0, 0, p)	
		else:
			self.X = x
			self.Z = z

	# @staticmethod
	def simplify_point(self):
		# x // z
		X = fp2_div(self.X, self.Z, self.p)
		Z = ComplexPoint(1, 0, self.p)

		return ECPoint(X,Z,self.p)

	def __eq__(self, other):
		if isinstance(other, ECPoint):
 			return self.X == other.X and self.Z == other.Z
		return False	

	def X(self):
		return self.X

	def Z(self):
		return self.Z

	def is_POIF(self):
		return self == ECPoint(x = ComplexPoint(0,0,self.p), p = self.p, y = ComplexPoint(1,0,self.p), z = ComplexPoint(0,0,self.p))

	def __str__(self):

		'''
		Converts a point P=(x,y) to P=(X,Z) 

		Params:
		P: point in affine coordinates

		Return:
		P in projectve coordinates
		'''
		
		return "(X:Z) = ({}:{})".format(self.X,self.Z)
