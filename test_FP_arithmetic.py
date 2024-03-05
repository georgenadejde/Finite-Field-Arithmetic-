import pytest
import unittest
from math import gcd
from FP_arithmetic import *
from hypothesis import given
from hypothesis.strategies import integers as rand


class TestUtilFunctions(unittest.TestCase):

    @given(rand(), rand())
    def gcd (self, a,b):
        assert Util.gcd(a, b) == math.gcd(b, a)

    @given(rand(), rand(), rand(2))
    def mod_power(self, x, p, m):
        result = Util.mod_power(x, p, m)
        assert result == pow(x,p,m)  # Compare with original python function

    @given(rand(), rand(2))
    def mul_inverse(self, other, mod):
        gcd = Util.gcd(other, mod)
        if gcd == 1:
            inverse = Util.mul_inverse(other, mod)
            # check the inverse is computed correctly
            assert (other * inverse) % mod == 1 
        else:
            result = Util.mul_inverse(other, mod)
            # check error message is printed
            assert result == "Inverse does not exist."

class TestFPArithmetic(unittest.TestCase):

    @given(rand(), rand(), rand(2))
    def test_finite_field_addition(self, a, b, mod):
        ModP = FiniteField(mod)
        a = ModP(a)
        b = ModP(b)
        result = a + b
        assert result.result == (a.elem + b.elem) % mod

    @given(rand(), rand(), rand(2))
    def test_finite_field_subtraction(self, a, b, mod):
        ModP = FiniteField(mod)
        a = ModP(a)
        b = ModP(b)
        result = a - b
        assert result.result == (a.elem - b.elem) % mod

    @given(rand(), rand(), rand(2))
    def test_finite_field_multiplication(self, a, b, mod):
        ModP = FiniteField(mod)
        a = ModP(a)
        b = ModP(b)
        result = a * b
        assert result.result == (a.elem * b.elem) % mod

    @given(rand(), rand(), rand(2))
    def test_finite_field_division(self, a, b, mod):
        ModP = FiniteField(mod)
        a = ModP(a)
        b = ModP(b)

        # Avoid division by zero and ensure division is possible
        if b.elem != 0 and Util.gcd(b.elem,mod) == 1:  
            result = a // b
            expected_result = (a.elem * Util.mul_inverse(b.elem, mod)) % mod
            assert result.result == expected_result

if __name__ == "__main__":
    unittest.main()