from Util import *
from ComplexPoint import *

def fp_add(x,y,p):
	'''
	Finite field addition

	Params:

	Return:
	Result finite field addition
	'''
	return (x + y) % p

def fp_dbl(x,p):

	return fp_add(x, x, p)

def fp_sub(x,y,p):
	'''
	Finite field subtraction

	Params:

	Return:
	finite field subtraction 
	'''

	return (x - y) % p

def fp_mul(x,y,p):
	'''
	Finite field multiplication

	Params:

	Return:
	Result finite field multiplication 
	'''
	return (x * y) % p

def fp_pow(x,pw,p):

	if pw == 0:
		return 1

	rez = fp_pow(x, pw // 2, p) % p
	square = fp_mul(rez,rez,p)
	if pw & 1:
		return (x * square) % p
	else:
		return square % p

def fp_inv(x,p):


	return Util.mul_inverse(x,p) % p

def fp_div(x,y,p):
	'''
	Finite field division

	Params:

	Return:

	'''
	return (x * fp_inv(y, p)) % p

def fp2_add(c1,c2,p):
	'''
	Finite field addition

	Params:
	c1: complex nr (a+bi)
	c2: complex nr (c+di)

	Return:
	Result finite field addition
	'''
	P = p
	return ComplexPoint(fp_add(c1.real, c2.real, P), fp_add(c1.imag, c2.imag, P), p)

def fp2_sub(c1,c2,p):
	'''
	Finite field subtraction

	Params:

	Return:
	finite field subtraction 
	'''
	
	return ComplexPoint(fp_sub(c1.real, c2.real, p), fp_sub(c1.imag, c2.imag, p), p)	

def fp2_mul(c1,c2,p):
	'''
	Finite field multiplication

	Params:

	Return:
	Result finite field multiplication 
	'''
	
	# (a+bi) * (c+di) = (ac - bd) + (ad + bc)i
	real_part = fp_sub(fp_mul(c1.real, c2.real, p), fp_mul(c1.imag, c2.imag, p), p)
	imag_part = fp_add(fp_mul(c1.real, c2.imag, p), fp_mul(c1.imag, c2.real, p), p) 

	return ComplexPoint(real_part, imag_part, p)

def fp2_pow(x,pw,p):

	if pw == 0:
		return 1

	rez = x
	while pw > 1:
		rez = fp2_mul(rez,x,p)
		pw = pw - 1

	return rez

def fp2_div(c1,c2,p):
	'''
	Finite field division

	Params:

	Return:

	'''
	
	# (a+bi) / (c+di) = (a+bi) (c-di) / (c^2 + d^2)


	c3 = ComplexPoint(c2.real, -c2.imag, p)
	N  = fp2_mul(c1,c3,p)
	
	r4 = fp_pow(c2.real,2,p)
	i4 = fp_pow(c2.imag,2,p) 
	R  = fp_add(r4,i4,p)
	D  = fp_inv(R,p)

	R  = fp2_mul(N,D,p)
	
	return  R


