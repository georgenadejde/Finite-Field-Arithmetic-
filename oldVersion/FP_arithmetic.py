from Util import *

def FiniteField(prime):
	class ModP:
		p = prime
		def __init__(self, elem):
			
			self.elem = elem % self.p 

		@property
		def result(self):
			return self.elem

		def __add__(self, other):
			'''
			Finite field addition

			Params:
			self: first point operand
			other: second point operand

			Return:
			Result finite field addition of self and other
			'''

			return ModP(self.elem + other.elem)

		def __sub__(self, other):
			'''
			Finite field subtraction

			Params:
			self: first point operand
			other: second point operand

			Return:
			Result finite field subtraction of self and other
			'''

			return ModP(self.elem - other.elem)

		def __mul__(self, other):
			'''
			Finite field multiplication

			Params:
			self: first point operand
			other: second point operand

			Return:
			Result finite field multiplication of self and other
			'''
			return ModP(self.elem * other.elem)
			

		def __floordiv__ (self, other):
			'''
			Finite field division

			Params:
			self: first point operand
			other: second point operand

			Return:
			self / other
			'''
			return ModP(self.elem * Util.mul_inverse(other.elem, ModP.p))

		
		def __str__(self):
			return str(self.elem) + " mod {}".format(ModP.p)

		

	return ModP
