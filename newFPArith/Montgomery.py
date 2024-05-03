# from Crypto.Util.number import inverse, getPrime
from FP_arith import * 
from Util import *
from ECPoint import *
import random
import cmath  


p = 1809251394333065553493296640760748573846642884164943950726130794095578185727 


# curve parameter
# A = ComplexPoint(751569296930863363598356131447593046593862768592437527947038742804342540969,
# 	            1798234580000369932369773894035327870837430306159012850472368713131657937742,
# 	            p)


xK = ComplexPoint(47740012980183545453926806418531833950371247694503085402799232463873975266
				  ,1609588601564889182030537691458495766012715929008094228092802417485560618033
				  , p)

yK = ComplexPoint(433448596440900625388421437356132198674684564930975145274932602841915121367,
	1448207449645360688773990653545547496727088082730664838864741114620683624430, p)


# point at infinity in projective coord. (X,Z)
POIF = ECPoint(ComplexPoint(1,0,p), p, z = ComplexPoint(0,0,p)) 

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

	if diff_P_Q.is_POIF() or diff_P_Q.is_T():
		return POIF
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

		return ECPoint(X, p, z=Z) 

def xDBL(P, a):
	'''
	Computes [2]P

	Params:
	P - point on the elliptic curve

	Return:
	[2]P
	'''	

	if P.is_POIF() or P.is_T():
		return POIF

	else:
		V_1  = fp2_add(P.X, P.Z, p)
		V_1  = fp2_mul(V_1, V_1, p)
		V_2  = fp2_sub(P.X, P.Z, p)
		V_2  = fp2_mul(V_2 ,V_2 ,p)
		x_2P = fp2_mul(V_1, V_2, p)

		V_1 =  fp2_sub(V_1, V_2, p)
		a2  =  ComplexPoint(2,0,p)
		c1  =  fp2_add(a,a2,p)
		c2  =  ComplexPoint(4, 0, p)

		V_3 = fp2_mul(fp2_div(c1, c2, p), V_1, p)
		V_3 = fp2_add(V_3, V_2, p)
		z_2P = fp2_mul(V_1, V_3, p)

		X = x_2P
		Z = z_2P

		return ECPoint(X, p, z= Z) 

def xDBLe(P, e, a):

	rez = P
	for i in range(e):
		rez = xDBL(rez, a)

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

def montgomery_ladder(k,P,a):
	'''
	Implementation of the Montgomery Ladder scalar
	multiplication algorithm.
	
	Params:
	P:  point in projective coordinates (X,Z)
	k:  scalar

	Return:
	[k]P in projective coordinates (X',Z')
	'''

	if P.is_POIF() or P.is_T():
		return POIF

	R_0 = P
	R_1 = xDBL(P,a)

	for bit in int_to_bin_list(k)[1:]:
		if bit == 0:

			R_1 = xADD(R_0, R_1, P) 
			R_0 = xDBL(R_0,a)

		else:
			R_0 = xADD(R_0,R_1, P)
			R_1 = xDBL(R_1,a)

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

def get_next_A(alpha,p):
	
	sq = fp2_mul(alpha,alpha,p) # alpha^2
	sq2 = fp2_add(sq,sq,p)   # 2*alpha^2
	c = ComplexPoint(1,0,p) # 1 + 0*i
	B = fp2_sub(c,sq2,p)    # 1-2*alpha^2
	R = fp2_add(B,B,p)      # 2*(1-2*alpha^2)

	return R

def iso2_curve(x,alpha,p):
 
	ax = fp2_mul(alpha,x,p)  # alpha*x
	c = ComplexPoint(1,0,p)  # 1 + 0*i
	n2 = fp2_sub(ax,c,p)     # alpha * x - 1 
	N = fp2_mul(x,n2,p)      # x * (alpha*x - 1)

	D = fp2_sub(x,alpha,p)   # x - alpha

	R = fp2_div(N,D,p)       # x * (alpha*x - 1) / (x-alpha)

	return R

def is_order_2(P,a):

	return xDBL(P,a).is_POIF()

def is_square(q):
	'''Legendre check for square
	
	[description]
	
	Args:
		q ([type]): [description]
	'''
	pw = fp_pow(q,(p-1) // 2, p)
	
	return (pw == 1) 

def get_norm(c):
	'''Compute N(c) = a^2 + b^2, where c = a+bi
	
	[description]
	
	Args:
		c ([type]): [description]
	'''

	a2 = fp_mul(c.real, c.real, p)
	b2 = fp_mul(c.imag, c.imag, p)

	norm = fp_add(a2,b2,p)

	return norm

def compute_y_squared(x,A):
	x2 = fp2_pow(x,2,p) 			# x^2 
	ax2 = fp2_mul(A,x2,p) 			# a*x^2
	x3 = fp2_pow(x,3,p) 			# x^3
	add1 = fp2_add(x3,ax2,p) 		# x^3 + a*x^2
	y_square = fp2_add(add1,x,p)    # x^3 + a*x^2 + x

	return y_square

def is_supersingular(A, num_trials=10):
    
    trial = 1
    while trial <= num_trials:
        # Get a random x coordinate on the EC
        xP = ComplexPoint(random.randint(2,p-1),random.randint(2,p-1),p)
        y_square = compute_y_squared(xP,A)
        norm = get_norm(y_square)
        # Check if x^3 + Ax^2 + x is a square  
        if is_square(norm):
            P = ECPoint(xP,p)
            ladder = montgomery_ladder(p*p-1,P,A)
            if ladder.is_POIF():
                trial = trial + 1
            else:
                return False
            
    return True

def isog_2e(a, S, e):
	
	while e:
		print(f"S.x = {S.X}")
		print(f"A = {a}")
		print(f"order = {e}")
		#check that a is supersingular
		assert is_supersingular(a)
		ker = xDBLe(S, e - 1, a) 
		print(f"ker = {ker}\n")
		print(f"XDBL = {xDBL(ker,a)}")

		ker.simplify_point()
		alpha = ker.X
		# verify order alpha is 2
		if is_order_2(ker,a):
			# get_next_A() needs alpha affine, e.g. simplified ker.X
			# instead, would want to do this projectively
			a = get_next_A(alpha, p)
			# print(f"j invariant = {j_invariant(a,p)}")
			e = e - 1
			if e  > 0:
				S.simplify_point()
				#this is now affine version
				S.X = iso2_curve(S.X, alpha, p)
				# iso2_curve_projective, with S projectve and alpha projective
				# ensure that (imX : imZ) simplifies to the same as (S.X : 1)
		else:
			print(f"\n{e}")
			print("Alpha is not of order 2.\n")
			break
		print("\n\n")
	return a


def main():

	'''
	QUESTIONS:
	- when checking if a point is of order 2, i get Z=0 but X is not 1. Is this still POIF?	
	- if i don't simplify the ker, the program does not work (yields a non-order 2 alpha for the next iter)
	'''


	A = ComplexPoint(751569296930863363598356131447593046593862768592437527947038742804342540969,
	            1798234580000369932369773894035327870837430306159012850472368713131657937742,
	            p)
	order = 52

	S = ECPoint(xK, p=p)

	rez = isog_2e(A,S,order)

	# print(S)
	# dblS = xDBLe(S, 51, A)
	# print(fp2_div(dblS.X, dblS.Z, p)) 	
	# print(dblS.simplify_point())

	print(f"My final A = {rez}\n")
	print(f"My j_invariant: {j_invariant(rez,p)}")

	x = ComplexPoint(733129935304444057399094841997377855449615348932625102326641103623966350086,
					622764407792741832064267776029604004876062748457690044116364384817066521173,p)
	z = ComplexPoint(971241277239417404462594040891046505615025646569885413196035817346782477360,
				1602963966652920992668623644809404193323285854490388641139684664564451659961,p)
	correctRez = fp2_div(x,z,p) 
	print(f"\n\nCORRECT j_invariant: {j_invariant(correctRez,p)}")

	print(f"\n\n ARE THEY THE SAME? {j_invariant(rez,p) == j_invariant(correctRez,p)}")

	sageA = ComplexPoint(1257067101235690800956486482453375172900760553591982339278706906628800031817, 
		390949035742497759430109018042411859860852776813963475965028075830395941209,p)

	print(f"SAGE result (j_invariant): {j_invariant(sageA,p)}")

main()