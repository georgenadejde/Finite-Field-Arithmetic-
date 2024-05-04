from Util import *

class Fp2Point:
    def __init__(self, real, imag, p=1):
        self.real = real
        self.imag = imag
        self.p = p

    def __eq__(self, other):
        if isinstance(other, Fp2Point):
            return self.real == other.real and self.imag == other.imag # and self.p == other.p

        return False

    def __str__(self):
        if self.imag == 0:
            return f"{self.real}"
        else:
            return f"{self.real} + i * {self.imag}"
