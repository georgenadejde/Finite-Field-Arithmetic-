from Util import *

def FiniteField(prime):
	class ModP:
		p = prime
		def __init__(self, elem):
			# print (elem, type(elem))
			# print(self.p, type(self.p))
			# print("\n\n")
			self.elem = elem % self.p 
			# self.elem = elem - (elem // self.p) * self.p 

		def __add__(self, other):
			'''
			Finite field addition

			Params:
			a: first point to add
			b: second point to add

			Return:
			Result finite field addition of a and b
			'''

			return ModP(self.elem + other.elem)

		def __sub__(self, other):
			'''
			Finite field subtraction

			Params:
			a: first point to add
			b: second point to add

			Return:
			Result finite field subtraction of a and b
			'''

			return ModP(self.elem - other.elem)

		def __mul__(self, other):
			'''
			Finite field multiplication

			Params:
			a: first point to add
			b: second point to add

			Return:
			Result finite field multiplication of a and b
			'''
			return ModP(self.elem * other.elem)
			

		def __floordiv__ (self, other):
			return ModP(self.elem * Util.mul_inverse(other.elem, ModP.p))

		def __str__(self):
			return str(self.elem) + " mod {}".format(ModP.p)

		@property
		def result(self):
			return self.elem

	return ModP


modP = FiniteField(7)

# print(modP(4) // modP(2))
# lol = modP(10)+modP(8)
# print(type(lol))
# print(type(modP))