from sympy import isprime

class Util:

	@staticmethod
	def isPrime(n):
		"""
		Check if a number is prime.

		n: The number to check
		
		Return: True if the number is prime, False otherwise
		"""
		if (n == 1):
			return False
		elif (n == 2):
			return True
		else:
			for x in range(2, n):
				if (n % x == 0):
					return False
				return True

		return False

	@staticmethod
	def gcd (a,b):
		"""
        Calculate the greatest common divisor (GCD) of two numbers.

        a: First number
        b: Second number
        
        Return: gcd(a, b)
        """

		while b != 0:
			r = a % b
			a = b
			b = r
		return a

	def mod_power(x, p, m):
		'''Calculate modular exponentiation (x^p % m).

        x: Base
        p: Exponent
        m: Modulus
        
        Return: Result of x^p % m
		'''

		if p == 0:
			return 1

		result = Util.mod_power(x, p // 2, m) % m
		square = (result * result) % m
		if p & 1:
			return (x * square) % m
		else:
			return square % m


	def mul_inverse(other, mod):
		"""
		Calculate the modular multiplicative inverse.

        other: Number for which inverse is to be calculated
        mod: Modulus
        
        Return: Modular multiplicative inverse of other modulo mod
		"""

		if(Util.gcd(other,mod) != 1 or not Util.isPrime(mod)):
			return "Inverse does not exist."
		else:
			return Util.mod_power(other, mod - 2, mod)
