import numpy as np 
from fractions import gcd
import random as rnd

class RationalNumber:
    def __init__(self, num, denom):
        gcdiv = gcd(num, denom)
        num = num / gcdiv
        denom = denom / gcdiv

        if denom < 0 and num > 0:
            num *= -1
            denom *= -1

        self._num = num
        self._denom = denom

    @staticmethod
    def Zero():
        return RationalNumber(0, 1)

    @staticmethod
    def One():
        return RationalNumber(1, 1)

    @staticmethod
    def create_random(low = -10, up = +10):
        num, denom = 0, 0

        while denom == 0:
            num, denom = rnd.randint(low, up), rnd.randint(low, up)

        return RationalNumber(num, denom)

    def Num(self):
        return self._num

    def Denom(self):
        return self._denom

    def __add__(self, o):
        if not isinstance(o, RationalNumber):
            raise TypeError("Expected rational number, got {0}".format(type(o)))
        
        return RationalNumber(self._num * o.Denom() + self._denom * o.Num(), self._denom * o.Denom())

    def __sub__(self, o):
        if not isinstance(o, RationalNumber):
            raise TypeError("Expected rational number, got {0}".format(type(o)))
        
        return RationalNumber(self._num * o.Denom() - self._denom * o.Num(), self._denom * o.Denom())

    def __mul__(self, o):
        if not isinstance(o, RationalNumber):
            raise TypeError("Expected rational number, got {0}".format(type(o)))
        
        return RationalNumber(self._num * o.Num(), self._denom * o.Denom())

    def __div__(self, o):
        if not isinstance(o, RationalNumber):
            raise TypeError("Expected rational number, got {0}".format(type(o)))
        
        return RationalNumber(self._num * o.Denom(), self._denom * o.Num())

    def __neg__(self):
        return RationalNumber(- self._num, self._denom)

    def __lt__(self, o):
        if not isinstance(o, RationalNumber):
            raise TypeError("Expected rational number, got {0}".format(type(o)))
        
        return self._num / self._denom < o.Num() / self.Denom()

    def __str__(self):
        if self._denom == 1:
            return str(self._num)
        elif self._num < 0:
            return "- " + str(- self._num) + "/" + str(self._denom)
        elif self._denom < 0:
            return "- " + str(self._num) + "/" + str(- self._denom)
        else:
            return str(self._num) + "/" + str(self._denom)

class RationalPolynomial:
    def __init__(self, coeff):
        coeff_red = [c for c in coeff if c is not RationalNumber.Zero()]
        self._coeff = coeff_red
        self._degree = len(self._coeff) - 1

    @staticmethod
    def create_random_rat_pol(degree):
        coeff = []

        for i in range(0, degree + 1):
            coeff.append(RationalNumber.create_random())

        return RationalPolynomial(coeff)

    def Coeff(self):
        return self._coeff

    def Degree(self):
        return self._degree

    def __add__(self, o):
        if not isinstance(self, RationalPolynomial):
            raise TypeError("Expected rational poylnomial, got {0}".format(type(o)))
        
        coeff_add = []
        coeff_longer = self._degree

        if o.Degree() > self._degree:
            coeff_longer = o.Coeff()

        for i, c in enumerate(coeff_longer):
            if i < len(coeff_longer) - self._degree:
                coeff_add.append(c)
            else:
                coeff_add.append(c + o.Coeff()[i - len(coeff_longer) + self._degree])

        return RationalPolynomial(coeff_add)

    def __sub__(self, o):
        if not isinstance(self, RationalPolynomial):
            raise TypeError("Expected rational poylnomial, got {0}".format(type(o)))
        
        return None

    def __mul__(self, o):
        if not isinstance(self, RationalNumber):
            raise TypeError("Expected rational poylnomial, got {0}".format(type(o)))
        
        return None

    def __div__(self, o):
        if not isinstance(self, RationalNumber):
            raise TypeError("Expected rational poylnomial, got {0}".format(type(o)))
        
        return None

    def __str__(self):
        pol_str = ""
        for i, c in enumerate(self._coeff):
            pol_str += "X^" + str(self._degree - i) + ": " + str(c) + "\n"

        return pol_str
    

class FiniteFieldExtension:
    def __init__(self, mini_pol):
        self._mini_pol = mini_pol

rand_pol_5 = RationalPolynomial.create_random_rat_pol(5)
rand_pol_8 = RationalPolynomial.create_random_rat_pol(8)
print(rand_pol_5 + rand_pol_8)




