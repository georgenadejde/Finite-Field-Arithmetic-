# from Crypto.Util.number import inverse, getPrime
from FP_arith import * 
from Util import *
from ECPoint import *
import random
import cmath  

# finite field modulo
p = 1809251394333065553493296640760748573846642884164943950726130794095578185727 

# point at infinity in projective coord. (X,Z)
POIF = ECPoint(Fp2Point(1,0), p, z = Fp2Point(0,0)) 

def int_to_bin_list(num):
	'''
	Converts an integer to its binary repr. in a list

	Params:
	num (int): int to be converted

	Return:
	binary representation of num in a list
	'''

	return [int(i) for i in [*bin(num).lstrip('0b')] ]

def xADD(P,Q, diff_P_Q):

	'''
	Adds two points found on a Montgomery curve

	Params:
	P (ECPoint): point on the Montgomery curve
	Q (ECPoint): point on the Montgomery curve
	diff_P_Q (ECPoint): P - Q

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

def xDBL(P, A):
	'''
	Computes [2]P

	Params:
	P (ECPoint): Point on the curve
	A (ECPoint): Curve parameter in proj. coord.

	Return:
	[2]P
	'''	

	if P.is_POIF() or P.is_T():
		return POIF

	else:
		# check if a is in projective coordinates
		# get x = X / Z
		if isinstance(A,ECPoint):
			A = A.normalize().X
		V_1  = fp2_add(P.X, P.Z, p)
		V_1  = fp2_mul(V_1, V_1, p)
		V_2  = fp2_sub(P.X, P.Z, p)
		V_2  = fp2_mul(V_2 ,V_2 ,p)
		x_2P = fp2_mul(V_1, V_2, p)

		V_1 =  fp2_sub(V_1, V_2, p)
		a2  =  Fp2Point(2,0)
		c1  =  fp2_add(A,a2,p)
		c2  =  Fp2Point(4, 0)

		V_3 = fp2_mul(fp2_div(c1, c2, p), V_1, p)
		V_3 = fp2_add(V_3, V_2, p)
		z_2P = fp2_mul(V_1, V_3, p)

		X = x_2P
		Z = z_2P

		return ECPoint(X, p, z= Z) 

def xDBLe(P, e, A):
	'''Computes [2^e]P
	
	Uses xDBL to scalar multiply P by 2^e
	
	Args:
		P (ECPoint): Point on the curve
		e (int): integer
		A (ECPoint): Curve parameter in proj. coord.
	
	Returns:
		[2^e]P
	'''
	rez = P
	for i in range(e):
		rez = xDBL(rez, A)

	return rez

def montgomery_ladder(k,P,A):
	'''
	Computes [k]P on curve E_A

	Uses the Montgomery Ladder scalar multiplication algorithm
	to compute [k]P.
	
	Params:
	k (int):  scalar
	P (ECPoint):  point in projective coordinates 
	A (ECPoint): curve parameter in proj. coordinates

	Return:
	[k]P in projective coordinates
	'''

	if P.is_POIF() or P.is_T():
		return POIF

	R_0 = P
	R_1 = xDBL(P,A)

	for bit in int_to_bin_list(k)[1:]:
		if bit == 0:

			R_1 = xADD(R_0, R_1, P) 
			R_0 = xDBL(R_0,A)

		else:
			R_0 = xADD(R_0,R_1, P)
			R_1 = xDBL(R_1,A)

	return R_0

def is_order_2(P, A):
	'''Check if P is a point of order 2 on curve defined
	by paramter A
	
	Computes xDBL(P) to verify if P has order 2 on Montgmery
	curve defined by parameter A.
	
	Args:
		P (ECPoint): Point on the curve
		A (ECPoint): Curve parameter
	
	Returns:
		True if P has order 2 on E_A, False otherwise.
	'''
	return xDBL(P, A).is_POIF()

def is_square(q):
	'''Check if an integer is a square in F_p
	
	Uses Legendre's symbol to check if q is a square
	by computing q^(p-1 / 2) % p
	
	Args:
		q ([int]): integer
	'''
	pw = fp_pow(q,(p-1) // 2, p)
	
	return (pw == 1) 

def get_norm(c):
	'''Compute N(c) = a^2 + b^2, where c = a+bi,
	   in finite field F_p
		
	Args:
		c (ComplexNumber): complex number
	'''

	a2 = fp_mul(c.real, c.real, p)   # a^2  
	b2 = fp_mul(c.imag, c.imag, p)   # b^2

	norm = fp_add(a2,b2,p)           # a^2 + b^2

	return norm

def compute_y_squared(x,A):

	'''Computes x^3 + Ax^2 + x 	
	
	
	Args:
		x (ECPoint): point on curve defined by curve param. A 
		A (ECPoint): Curve parameter in proj. coord.
	
	Returns:
		y^2 = x^3 + Ax^2 + x
	'''

	# check if a is in projective coordinates
	# get x = X / Z
	if isinstance(A,ECPoint):
		A = A.normalize().X
	
	x2 = fp2_pow(x,2,p) 			# x^2 
	ax2 = fp2_mul(A,x2,p) 			# a*x^2
	x3 = fp2_pow(x,3,p) 			# x^3
	add1 = fp2_add(x3,ax2,p) 		# x^3 + a*x^2
	y_square = fp2_add(add1,x,p)    # x^3 + a*x^2 + x

	return y_square

def is_supersingular(A, num_trials=10):
    '''Check if a curve is supersingular
    
    Runs a check on a Montgomery curve with parameter A
    on whether the curve is supersingular by computing
    [p^2-1]P on random points P on the curve.
    
    Args:
    	A (ECPoint): curve parameter in proj. coordinates
    	num_trials (number): number of trials (default: `10`)
    
    Returns:
    	True if the curve is supersingular, False otherwise.
    '''
    trial = 1
    while trial <= num_trials:
        # Get a random x coordinate on the EC
        xP = Fp2Point(random.randint(2,p-1),random.randint(2,p-1))
        # Computes x^3 + Ax^2 + x
        y_square = compute_y_squared(xP,A)
        # Computes the norm of y^2
        norm = get_norm(y_square)
        # Check if x^3 + Ax^2 + x is a square  
        if is_square(norm):
        	# Create the random point
            P = ECPoint(xP,p)
            # Computes [p^2-1]P
            ladder = montgomery_ladder(p*p-1,P,A)
            if ladder.is_POIF():
                trial = trial + 1
            else:
                return False
            
    return True

def j_invariant(A):
	'''
	Computes the j-invariant of a Montgomery curve
	
	Params:
	a: Montgomery coefficient in projective coordinates

	Return:
	j(E_a)
	'''

	# set the finite field modulo (for testing purposes)
	# can also add p as a parameter
	p = A.p

	c256 = Fp2Point(256,0)    # 256 + 0*i
	c3 = Fp2Point(3,0)        # 3 + 0*i
	c4 = Fp2Point(4,0)        # 4 + 0*i

	ax2 = fp2_mul(A.X,A.X,p)        # A.X^2
	az2 = fp2_mul(A.Z,A.Z,p)        # A.Z^2
	az4 = fp2_mul(az2,az2,p)        # A.Z^4
	
	az2_3 = fp2_mul(c3, az2, p)		# 3*A.Z^2
	az2_4 = fp2_add(az2_3, az2, p)  # 4*A.Z^2
	
	subN = fp2_sub(ax2,az2_3,p)     # A.X^2 - 3*A.Z^2
	subD = fp2_sub(ax2,az2_4,p)     # A.X^2 - 4*A.Z^2

	subN3 = fp2_pow(subN,3,p)       # (A.X^2 - 4*A.Z^2)^3
	X = fp2_mul(c256, subN3, p)		# 256 * (A.X^2 - 4*A.Z^2)^3

	Z = fp2_mul(az4, subD, p)       # A.Z^4 * (A.X^2 - 4A.Z^2)

	# get affine coordinates
	# can also return (X:Z) for projective coord.
	R = fp2_div(X,Z,p)              

	return R 

def get_next_A(alpha):
	
	'''Computes the A parameter of the codomain's Mont. curve
	
	Uses Craig's formula to compute the A parameter of the 
	isogeny's codomain Montgomery curve
	
	Args:
		alpha (ECPoint): kernel in proj. coord.
	
	Returns:
		A parameter in proj. coordinates
	'''

	# set the finite field modulo (for testing purposes)
	# can also add p as a parameter
	p = alpha.p
	
	sqX = fp2_mul(alpha.X,alpha.X,p)    # alpha.X^2
	sqZ = fp2_mul(alpha.Z,alpha.Z,p)    # alpha.Z^2
	
	sqX2 = fp2_dbl(sqX,p)           # 2*alpha.X^2
	sub =  fp2_sub(sqZ,sqX2,p)          # alpha.Z^2 - 2*alpha.X^2
	
	X = fp2_dbl(sub,p)              # 2(alpha.Z^2 - 2*alpha.X^2)

	Z = sqZ

	# returns A in projective coordinates
	R = ECPoint(X,p,z=Z)

	return R

def iso2_curve(x,alpha):
	'''Computes a 2-isogeny
	
	Uses a map to compute a 2-isogeny using a point x on the curve and 
	the kernel alpha
	
	Args:
		x (ECPoint): a point of the curve in the domain of the isogeny
		alpha (ECPoint): kernel of the curve in the domain of the isogeny
		                 in projective coordinates.
	
	Returns:
		The mapped point of x on the curve in the codomain of the isogeny
	'''
	SxAx = fp2_mul(x.X, alpha.X, p)   # S.X * Alpha.X 
	SzAz = fp2_mul(x.Z, alpha.Z, p)   # S.Z * Alpha.Z
	SxAz = fp2_mul(x.X, alpha.Z, p)   # S.X * Alpha.Z
	SzAx = fp2_mul(x.Z, alpha.X, p)   # S.Z * Alpha.X

	P = fp2_sub(SxAx,SzAz,p)          # S.X * Alpha.X - S.Z * Alpha.Z          
	Q = fp2_sub(SxAz,SzAx,p)          # S.X * Alpha.Z - S.Z * Alpha.X

	X = fp2_mul(x.X, P, p)            # S.X * (S.X * Alpha.X - S.Z * Alpha.Z)
	Z = fp2_mul(x.Z, Q, p)            # S.Z * (S.X * Alpha.Z - S.Z * Alpha.X)

	# return in projective coordinates (X:Z)
	R = ECPoint(X,p,z=Z)

	return R

def isog_2e(A, S, e):
	'''Computes a 2^e isogeny
	
	Uses a chain of e 2-isogenies to compute a 2^e isogeny and returns the A parameter of the destination curve 		
	
	Args:
		a (ECPoint): Montgomery curve parameter 
		             in projective coordinates 	
		S (ECPoint): A point on curve defined by a, of order 2^e 
		             in projective coordinates.
		e (int): order of the isogeny
	
	Returns:
		Parameter a of the destination curve in proj. coord. (X:1)
	'''
	while e:
		
		#check that a is supersingular
		assert is_supersingular(A)
		
		ker = xDBLe(S, e - 1, A) 

		# verify order ker is 2
		if is_order_2(ker,A):
			
			# gets the isogeny codomain's a parameter
			A = get_next_A(ker)
			
			e = e - 1
			if e  > 0:
				# gets the next point of order 2^e-1
				S = iso2_curve(S, ker)
		else:
			print("Kernel is not of order 2.\n")
			break

	# check we arrived on a supersingular curve	
	assert is_supersingular(A)

	# return A in simplified (affine) coordinates
	return A.normalize()

def main():

	A = ECPoint(Fp2Point(751569296930863363598356131447593046593862768592437527947038742804342540969,
	            1798234580000369932369773894035327870837430306159012850472368713131657937742),p)
	
	xK = Fp2Point(47740012980183545453926806418531833950371247694503085402799232463873975266
				  ,1609588601564889182030537691458495766012715929008094228092802417485560618033)

	order = 52

	S = ECPoint(xK, p)

	rez = isog_2e(A,S,order)


	print(f"My final A = {rez}\n")
	print(f"My j_invariant: {j_invariant(rez)}")

	x = Fp2Point(733129935304444057399094841997377855449615348932625102326641103623966350086,
					622764407792741832064267776029604004876062748457690044116364384817066521173)
	z = Fp2Point(971241277239417404462594040891046505615025646569885413196035817346782477360,
				1602963966652920992668623644809404193323285854490388641139684664564451659961)
	
	correctRez = ECPoint(x,p,z=z)
	print(f"\n\nCORRECT j_invariant: {j_invariant(correctRez)}")

	print(f"\n\n ARE THEY THE SAME? {j_invariant(rez) == j_invariant(correctRez)}")

main()