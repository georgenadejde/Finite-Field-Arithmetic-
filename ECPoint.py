from Util import *


class ECPoint:
	def __init__(self, x, y = 0, z = 1, modP = None):
		
		# y is not needed for montgomery curves
		self.y = None

		# poif in projective coordinates
		if not x and not z:
			self.x = 1
			self.z = z	
		else:
			self.x = x
			self.z = z

		# added just in case we want to link a point to a field
		# and automatically simplify the point
		if modP is not None:
			self.simplify_point(self, modP)

	# @staticmethod
	def simplify_point(self, modP):
		result = modP(self.x) // modP(self.z)
		self.x = result.elem
		self.z = 1
		return self
		# return ECPoint(result.elem, z = 1)

	def __eq__(self, other):
		if isinstance(other, ECPoint):
 			return self.x == other.x and self.z == other.z
		return False	

	def x(self):
		return self.x

	def z(self):
		return self.z

	def is_POIF(self):
		return self == ECPoint(0,1,0)

	def __str__(self):

		'''
		Converts a point P=(x,y) to P=(X,Z) 

		Params:
		P: point in affine coordinates

		Return:
		P in projectve coordinates
		'''
		
		return "({}:{})".format(self.x,self.z)
