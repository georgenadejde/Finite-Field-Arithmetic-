from Util import *
from ComplexPoint import *

def fp_add(x,y,p):
	'''Finite field F_p addition
	
	Args:
		x (int): integer
		y (int): integer
		p (int): modulo
	
	Returns:
		(x+y)%p
	'''
	return (x + y) % p

def fp_dbl(x,p):
	'''Finite field F_p doubling
	
	Args:
		x (int): integer
		p (int): modulo
	
	Returns:
		2*x % p
	'''

	return fp_add(x, x, p)

def fp_sub(x,y,p):
	'''Finite field F_p subtraction
	
	Args:
		x (int): integer
		y (int): integer
		p (int): modulo
	
	Returns:
		(x-y)%p
	'''

	return (x - y) % p

def fp_mul(x,y,p):
	'''Finite field F_p multiplication
	
	Args:
		x (int): integer
		y (int): integer
		p (int): modulo
	
	Returns:
		(x+y)%p
	'''

	return (x * y) % p

def fp_pow(x,pw,p):
	'''Finite field F_p pow function
	
	Args:
		x (int): integer
		pw (int): integer
		p (int): modulo
	
	Returns:
		x^pw % p
	'''
	if pw == 0:
		return 1

	rez = fp_pow(x, pw // 2, p) % p
	square = fp_mul(rez,rez,p)
	if pw & 1:
		return (x * square) % p
	else:
		return square % p

def fp_inv(x,p):

	'''Computes the multiplicative inverse in F_p
	
	Args:
		x (int): integer
		p (int): modulo
	
	Returns:
		y so that (x * y) % p = 1
	'''

	return Util.mul_inverse(x,p) % p

def fp_div(x,y,p):
	'''Finite field F_p division
	
	Args:
		x (int): integer
		y (int): integer
		p (int): modulo
	
	Returns:
		(x / y) % p
	'''
	return (x * fp_inv(y, p)) % p

def fp2_add(c1,c2,p):
	'''
	Finite field extension F_p^2 addition

	Params:
	c1: complex nr (a+bi)
	c2: complex nr (c+di)
	p: modulo

	Return:
		(c1 + c2) % p
	'''

	return ComplexPoint(fp_add(c1.real, c2.real, p), fp_add(c1.imag, c2.imag, p))

def fp2_dbl(x,p):
	'''Finite field extension F_p2 doubling
	
	Args:
		x (int): complex number a+bi
		p (int): modulo
	
	Returns:
		2*x % p
	'''

	return fp2_add(x, x, p)

def fp2_sub(c1,c2,p):
	'''
	Finite field extension F_p^2 subtraction

	Params:
	c1: complex nr (a+bi)
	c2: complex nr (c+di)
	p: modulo

	Return:
		(c1 - c2) % p
	'''
	
	return ComplexPoint(fp_sub(c1.real, c2.real, p), fp_sub(c1.imag, c2.imag, p))	

def fp2_mul(c1,c2,p):
	'''
	Finite field extension F_p^2 multiplication

	Params:
	c1: complex nr (a+bi)
	c2: complex nr (c+di)
	p: modulo

	Return:
		(c1 * c2) % p
	'''
	
	# (a+bi) * (c+di) = (ac - bd) + (ad + bc)i
	real_part = fp_sub(fp_mul(c1.real, c2.real, p), fp_mul(c1.imag, c2.imag, p), p)
	imag_part = fp_add(fp_mul(c1.real, c2.imag, p), fp_mul(c1.imag, c2.real, p), p) 

	return ComplexPoint(real_part, imag_part)

def fp2_pow(x,pw,p):

	'''
	Finite field extension F_p^2 pow function

	Params:
	x: complex nr (a+bi)
	pw: exponent
	p: modulo

	Return:
		(x^pw) % p
	'''

	if pw == 0:
		return 1

	rez = x
	while pw > 1:
		rez = fp2_mul(rez,x,p)
		pw = pw - 1

	return rez

def fp2_div(c1,c2,p):
	'''
	Finite field extension F_p^2 division

	Params:
	c1: complex nr (a+bi)
	c2: complex nr (c+di)
	p: modulo

	Return:
		(c1 / c2) % p
	'''
	
	# (a+bi) / (c+di) = (a+bi) (c-di) / (c^2 + d^2)

	c3 = ComplexPoint(c2.real, -c2.imag)
	N  = fp2_mul(c1,c3,p)
	
	r4 = fp_pow(c2.real,2,p)
	i4 = fp_pow(c2.imag,2,p) 
	R  = fp_add(r4,i4,p)
	D  = fp_inv(R,p)

	R  = fp2_mul(N,D,p)
	
	return  R


