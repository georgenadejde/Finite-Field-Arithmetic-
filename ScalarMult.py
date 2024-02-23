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
	
	# print("LOL\n")
	# print(type(modP(P.x)))
	# print(type(modP(P.z)))
	# print("\n\n")

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

	# print((z))
	# print((x))
	# print("\n\n")

	return ECPoint(X, 0, Z) 

def xDBL(P):
	
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
	Z = z_2P.result

	return ECPoint(X, 0, Z) 


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
	
	for bit in int_to_bin_list(k):
		# print("LOL\n")
		print("R_0=",R_0)
		print("R_1=",R_1)
		# print(type(P))
		

		if bit == 0:
			print("bit: 0\n")
			R_1 = xADD(R_0, R_1, P) 
			R_0 = xDBL(R_0)
		else:
			print("bit: 1\n")
			R_0 = xADD(R_0,R_1, P)
			R_1 = xDBL(R_1)

		print("\n\n")
	
	print("R_0=",R_0)
	print("R_1=",R_1)	

	return R_0
		

def main():
	print (int_to_bin_list(17))

	# xADD(ECPoint(1,0,4), ECPoint(10,0,7), ECPoint(1,0,4))
	# xDBL(ECPoint(1,0,4))
	P = ECPoint(387)
	A = 4
	rez = montgomery_ladder(2,P)
	print("Results")
	print(xDBL(ECPoint(387)))
	print("Rez:",rez)
	# print(modP(rez.x) // modP(rez.z) )



main()