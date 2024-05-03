import pytest
import unittest
from hypothesis import given
from hypothesis.strategies import integers as rand
from FP_arith import *
from ComplexPoint import *
from ECPoint import *
from math import gcd
from sympy import isprime
from Montgomery import *

p = 431

class TestUtilFunctions(unittest.TestCase):

    @given(rand(2), rand(2))
    def test_gcd(self, a, b):
        assert Util.gcd(a, b) == gcd(b, a)

    @given(rand(2), rand(2), rand(2))
    def test_mod_power(self, x, p, m):
        result = (Util.mod_power(x, p, m)) % m
        assert result == pow(x,p,m)  # Compare with original python function

    @given(rand(2), rand(2))
    def test_mul_inverse(self, other, mod):
        if isprime(mod) and Util.gcd(other, mod) == 1:
            inverse = Util.mul_inverse(other, mod)
            # check the inverse is computed correctly
            assert (other * inverse) % mod == 1 

class TestFPArithmetic(unittest.TestCase):

    @given(rand(), rand(), rand(2))
    def test_finite_field_addition(self, a, b, mod):
        result = fp_add(a, b, mod)
        assert result == (a + b) % mod

    @given(rand(), rand(2))
    def test_finite_field_doubling(self, x, mod):
        result = fp_dbl(x, mod)
        assert result == (x + x) % mod

    @given(rand(), rand(), rand(2))
    def test_finite_field_subtraction(self, a, b, mod):
        result = fp_sub(a, b, mod)
        assert result == (a - b) % mod

    @given(rand(), rand(), rand(2))
    def test_finite_field_multiplication(self, a, b, mod):
        result = fp_mul(a, b, mod)
        assert result == (a * b) % mod

    @given(rand(min_value=0), rand(min_value=0), rand(min_value=2))
    def test_finite_field_power(self, x, pw, mod):
        result = fp_pow(x, pw, mod)
        assert result == pow(x, pw, mod)

    @given(rand(), rand(2))
    def test_finite_field_inverse(self, x, mod):
        gcd_val = Util.gcd(x, mod)
        if isprime(mod) and gcd_val == 1:
            result = fp_inv(x, mod)
            assert (x * result) % mod == 1
    
    @given(rand(), rand(), rand(2))
    def test_finite_field_division(self, a, b, mod):
        gcd_val = Util.gcd(b, mod)
        if isprime(mod) and gcd_val == 1:
            result = fp_div(a, b, mod)
            assert result == (a * Util.mul_inverse(b, mod)) % mod
        
    @given(rand(), rand(), rand(), rand(), rand(2))
    def test_fp2_addition(self, a, b, c, d, mod):
        result = fp2_add(ComplexPoint(a, b), ComplexPoint(c, d), mod)
        assert result == ComplexPoint((a + c) % mod, (b + d) % mod)

    @given(rand(), rand(), rand(), rand(), rand(2))
    def test_fp2_subtraction(self, a, b, c, d, mod):
        result = fp2_sub(ComplexPoint(a, b), ComplexPoint(c, d), mod)
        assert result == ComplexPoint((a - c) % mod, (b - d) % mod)

    @given(rand(), rand(), rand(), rand(), rand(2))
    def test_fp2_multiplication(self, a, b, c, d, mod):
        result = fp2_mul(ComplexPoint(a, b), ComplexPoint(c, d), mod)
        real_part = (a * c - b * d) % mod
        imag_part = (a * d + b * c) % mod
        assert result == ComplexPoint(real_part, imag_part)

    # @given(rand(min_value=0, max_value=10), rand(min_value=0, max_value=10), rand(min_value=0, max_value=10), rand(min_value=2))
    # def test_fp2_power(self, a, b, pw, mod):
    #     c  = ComplexPoint(a,b,mod)
    #     rez = C
    #     if 
    #     result = fp2_pow(c, pw, mod)
    #     assert result == pow(x, pw, mod)

    @given(rand(), rand(), rand(), rand(), rand(2))
    def test_fp2_division(self, a, b, c, d, mod):
        c1 = ComplexPoint(a, b)
        c2 = ComplexPoint(c, d)
        gcd_val = Util.gcd(c2.real, mod)
        if isprime(mod*mod) and gcd_val == 1:
            result = fp2_div(c1, c2, mod)
            assert result == ComplexPoint((a * c2.real + b * c2.imag) % mod, (b * c2.real - a * c2.imag) % mod)

class TestIsogenyFormulas(unittest.TestCase):
    
    def test_dbl(self):
        pass
        # pg. 15

        # A = ComplexPoint(423,329)
        # print(xxDBLe(ECPoint(ComplexPoint(79,271),p), 3,A))
        # assert xxDBLe(ECPoint(ComplexPoint(79,271),p), 3,A) == ECPoint(ComplexPoint(37,18), p)

        # A = ComplexPoint(132,275)
        # # pg. 15
        # assert xxDBLe(ECPoint(ComplexPoint(111,36), p), 2) == ECPoint(ComplexPoint(49,7), p)

        # A = ComplexPoint(76,273)
        # # pg. 15
        # print(xxDBLe(ECPoint(ComplexPoint(374,274), p), 1,A))
        # assert xxDBLe(ECPoint(ComplexPoint(374,274), p), 1,A) == ECPoint(ComplexPoint(27,245), p)


    def test_j_invariant(self):
        # pg. 7 Craig
        print(j_invariant(ECPoint(ComplexPoint(415,0),p)))
        assert j_invariant(ECPoint(ComplexPoint(415,0),p)) == ComplexPoint(189,0)
        
        # pg. 6 Craig
        assert j_invariant(ECPoint(ComplexPoint(423,102),p)) == ComplexPoint(190,344)

        # pg. 6 Craig
        assert j_invariant(ECPoint(ComplexPoint(161,208),p)) == ComplexPoint(304,364)

        # pg.3 Craig
        assert j_invariant(ECPoint(ComplexPoint(162,172),p)) == ComplexPoint(304,364)

        # pg. 14 Craig
        assert j_invariant(ECPoint(ComplexPoint(423,329),p)) == ComplexPoint(190,87)

    def test_compute_A(self):

        # pg. 6
        alpha = ECPoint(ComplexPoint(68,350),p)
        newA = get_next_A(alpha)
        assert j_invariant(newA) == ComplexPoint(190,344)

        # pg. 15
        alpha = ECPoint(ComplexPoint(37,18),p)
        newA = get_next_A(alpha)
        assert j_invariant(newA) == ComplexPoint(107,0)

        # pg. 15
        alpha = ECPoint(ComplexPoint(49,7),p)
        newA = get_next_A(alpha)
        assert j_invariant(newA) == ComplexPoint(190, 344)

    

if __name__ == "__main__":
    unittest.main()
