# from Crypto.Util.number import inverse, getPrime
from FP_arith import * 
from Util import *
from ECPoint import *


p = 1809251394333065553493296640760748573846642884164943950726130794095578185727 


# curve parameter
A = ComplexPoint(751569296930863363598356131447593046593862768592437527947038742804342540969,
	            1798234580000369932369773894035327870837430306159012850472368713131657937742,
	            p)


xK = ComplexPoint(47740012980183545453926806418531833950371247694503085402799232463873975266
				  ,1609588601564889182030537691458495766012715929008094228092802417485560618033
				  , p)

# point at infinity in projective coord. (X,Z)
POIF = ECPoint(x = ComplexPoint(0,0,p), p = p, y = ComplexPoint(1,0,p), z = ComplexPoint(0,0,p)) 

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

	if diff_P_Q.is_POIF():
		return ECPoint(x=ComplexPoint(0,0,p), p=self.p, z=ComplexPoint(0,0,p)) 
	else:
		V_0 = fp2_add(P.X, P.Z, p)
		V_1 = fp2_sub(Q.X, Q.Z, p)
		V_1 = fp2_mul(V_1, V_0, p)
		V_0 = fp2_sub(P.X, P.Z, p)
		V_2 = fp2_add(Q.X, Q.Z, p)
		V_2 = fp2_mul(V_2, V_0, p)
		V_3 = fp2_add(V_1, V_2, p)
		V_3 = fp2_mul(V_3, V_3, p)
		V_4 = fp2_sub(V_1, V_2, p)
		V_4 = fp2_mul(V_4, V_4, p)

		x_rez = fp2_mul(diff_P_Q.Z, V_3, p)
		z_rez = fp2_mul(diff_P_Q.X, V_4, p)

		X = x_rez
		Z = z_rez

		return ECPoint(X, None, Z) 

def xDBL(P):
	'''
	Computes [2]P

	Params:
	P - point on the elliptic curve

	Return:
	[2]P
	'''
	
	V_1 = fp2_add(P.X, P.Z, p)
	V_1 = fp2_pow(V_1, 2, p)
	V_2 = fp2_sub(P.X, P.Z, p)
	V_2 = fp2_pow(V_2 ,2 ,p)
	x_2P = fp2_mul(V_1, V_2, p)

	V_1 = fp2_sub(V_1, V_2, p)

	a2 =  ComplexPoint(2,0,p)
	c1 =  fp2_add(A,a2,p)
	c2 =  ComplexPoint(4, 0, p)

	V_3 = fp2_mul(fp2_div(c1, c2, p), V_1, p)
	V_3 = fp2_add(V_3, V_2, p)
	z_2P = fp2_mul(V_1, V_3, p)

	X = x_2P

	if P.is_POIF():
		Z = 0
	else:
		Z = z_2P

	return ECPoint(X, p, y=None, z= Z) 

def xDBLe(P,e):

	rez = P
	for i in range(0, e):
		rez = xDBL(rez)

	return rez

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

def j_invariant(a,p):
	'''
	Computes the j-invariant of a Montgomery curve
	
	Params:
	a: Montgomery coefficient
	p: modulo of the finite field

	Return:
	j(E_a)
	'''

	V = fp2_mul(a,a,p) 			# a^2
	B = fp2_sub(V,3,p) 			# a^2 - 3
	C = fp2_pow(B,3,p) 			# (a^2-3)^3
	c1 = ComplexPoint(256,0,p)
	N = fp2_mul(c1, C, p)		# 256 * (a^2-3)^3

	c2 = ComplexPoint(4,0,p)
	D = fp2_sub(V,c2,p)         # a^2-4
	
	R = fp2_div(N,D,p)		    # [ 256 * (a^2-3)^3 ] / (a^2-4)

	return R

def get_A(alpha,p):
	
	sq = fp2_pow(alpha,2,p) # alpha^2
	sq2 = fp2_mul(sq,2,p)   # 2*alpha^2
	c = ComplexPoint(1,0,p) # 1 + 0*i
	B = fp2_sub(c,sq2,p)    # 1-2*alpha^2
	R = fp2_mul(2,B,p)      # 2*(1-2*alpha^2)

	return R

def iso2_curve(x,alpha,p):
 
	ax = fp2_mul(alpha,x,p)  # alpha*x
	c = ComplexPoint(-1,0,p) # -1 + 0*i
	n2 = fp2_sub(ax,c,p)     # alpha * x - 1 
	N = fp2_mul(x,n2,p)      # x * (alpha*x - 1)

	D = fp2_sub(x,alpha,p)   # x - alpha

	R = fp2_div(N,D,p)       # x * (alpha*x - 1) / (x-alpha)

	return R

def iso2_eval(x):
	pass

def isog_2e(a, S, e):
	
	while e:
		ker = xDBLe(S, e - 1) 
		alpha = ker.X
		a = get_A(alpha, p)
		e = e - 1

	return a


def main():

	'''
	QUESTIONS:

	'''
	A = ComplexPoint(751569296930863363598356131447593046593862768592437527947038742804342540969,
	            1798234580000369932369773894035327870837430306159012850472368713131657937742,
	            p)
	order = 52

	S = ECPoint(xK, p=p)
	

	A = isog_2e(A,S,order)
	
	print(j_invariant(A,p))


	c = ComplexPoint(733129935304444057399094841997377855449615348932625102326641103623966350086,
					622764407792741832064267776029604004876062748457690044116364384817066521173,p)

	print(f"\n\nCORRECT: {j_invariant(c,p)}")

	print(f"\n\n ARE THEY THE SAME? {j_invariant(A,p) == j_invariant(c,p)}")

	# P = ECPoint(ComplexPoint(387,0,p),p)
	# lad = montgomery_ladder(4,P)

	# print(lad.simplify_point())
	# print(modP(lad.x) // modP(lad.z))


	# rez = xADD(P,POIF,P)



main()