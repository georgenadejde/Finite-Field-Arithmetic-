# from Crypto.Util.number import inverse, getPrime
from FP_arithmetic import * 
from Util import *
from ECPoint import *

# curve parameters
A = 4
B = 1

p = 419 
# p = getPrime(37) # return a random 37 bit prime 
POIF = ECPoint(0,1,0)	# point at infinity in projective coord. (X,Z)

# Define an alias function for Util.mul_inverse
Inverse = Util.mul_inverse

# Define the modulo
modP = FiniteField(p)
	
def xADD(P,Q, diff_P_Q):

	'''
	Adds two elliptic curve points from a Montgomery curve

	Params:
	P: point on the elliptic curve
	Q: point on the elliptic curve
	diff_P_Q: P - Q

	Return:
	P + Q
	'''

	if diff_P_Q.is_POIF:
		return ECPoint(0, None, 0) 
	else:
		V_0 = modP(P.x) + modP(P.z)
		V_1 = modP(Q.x) - modP(Q.z)
		V_1 = V_1 * V_0
		V_0 = modP(P.x) - modP(P.z)
		V_2 = modP(Q.x) + modP(Q.z)
		V_2 = V_2 * V_0
		V_3 = V_1 + V_2
		V_3 = V_3 * V_3
		V_4 = V_1 - V_2
		V_4 = V_4 * V_4

		x_rez = (modP(diff_P_Q.z) * V_3)
		z_rez = (modP(diff_P_Q.x) * V_4)

		X = x_rez.result
		Z = z_rez.result

		return ECPoint(X, None, Z) 

def xDBL(P):
	'''
	Computes [2]P

	Params:
	P - point on the elliptic curve

	Return:
	[2]P
	'''
	
	V_1 = modP(P.x) + modP(P.z)
	V_1 = V_1 * V_1
	V_2 = modP(P.x) - modP(P.z)
	V_2 = V_2 * V_2
	x_2P = V_1 * V_2

	V_1 = V_1 - V_2
	V_3 = (modP(A+2) // modP(4)) * V_1
	V_3 = V_3 + V_2
	z_2P = V_1 * V_3

	X = x_2P.result

	if P.is_POIF():
		Z = 0
	else:
		Z = z_2P.result	

	return ECPoint(X, None, Z) 


def int_to_bin_list(num):
	'''
	Converts an integer to its binary repr. in a list

	Params:
	num: int to be converted

	Return:
	binary representation of num in a list
	'''

	return [int(i) for i in [*bin(num).lstrip('0b')] ]

def montgomery_ladder(k,P):
	'''
	Implementation of the Montgomery Ladder scalar
	multiplication algorithm.
	
	Params:
	P:  point in projective coordinates (X,Z)
	k:  scalar

	Return:
	[k]P in projective coordinates (X',Z')
	'''

	R_0 = P
	R_1 = xDBL(P)
	
	for bit in int_to_bin_list(k)[1:]:

		if bit == 0:
			R_1 = xADD(R_0, R_1, P) 
			R_0 = xDBL(R_0)
		else:
			R_0 = xADD(R_0,R_1, P)
			R_1 = xDBL(R_1)


	return R_0
		

def main():

	'''
	QUESTIONS:

	- What is T precisely?
	- [0]P ?
	- since we only work with the x-coordinate,
	  why is the ladder not returning the final
	  value of it? now, it misses the step where
	  it should be multiplied by the inverse of z
	- added a simplify_point method, but not sure if needed,
	  because the algo in the paper does not simplify the point
	'''

	P = ECPoint(387)
	lad = montgomery_ladder(4,P)

	print(lad.simplify_point(modP))
	print(modP(lad.x) // modP(lad.z))

	rez = xADD(P,POIF,P)
	# print(rez.x // rez.z)


main()