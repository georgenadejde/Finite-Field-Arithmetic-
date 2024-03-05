class Util:

	@staticmethod
	def gcd (a,b):
		"""
        Calculate the greatest common divisor (GCD) of two numbers.

        :param a: First number
        :param b: Second number
        :return: GCD of a and b
        """
		while b != 0:
			r = a % b
			a = b
			b = r
		return a

	def mod_power(x, p, m):
		'''Calculate modular exponentiation (x^p % m).

        :param x: Base
        :param p: Exponent
        :param m: Modulus
        :return: Result of x^p % m
		'''
		if p == 0:
			return 1

		result = Util.mod_power(x, p // 2, m) % m
		square = (result * result) % m
		if p & 1:
			return (x * square)
		else:
			return square


	def mul_inverse(other, mod):
		"""
		 Calculate the modular multiplicative inverse.

        :param other: Number for which inverse is to be calculated
        :param mod: Modulus
        :return: Modular multiplicative inverse of other modulo mod
		"""

		if(Util.gcd(other,mod) != 1):
			return "Inverse does not exist."
		else:
			return Util.mod_power(other, mod - 2, mod)
	