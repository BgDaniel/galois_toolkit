import numpy.polynomial.polynomial as poly
from numpy.polynomial.polynomial import Polynomial as P
from sympy.combinatorics.named_groups import SymmetricGroup
import mpmath
from mpmath import *
from random import randint

from GaloisResolvent import *
from Sym5 import Sym5
from Sym4 import Sym4
from Sym3 import Sym3


MAX_ITERATIONS = 5000
PRECISION = 0.00000001

class PolyEquation:

    @property
    def Coefficients(self):
        return self._coefficients

    @property
    def Degree(self):
        return self._degree

    @property
    def Roots(self):
        return self._roots

    def __init__(self, coefficients, dps = 50):
        self._coefficients = coefficients
        self._polynom = P(coefficients) 
        self._degree = len(coefficients) - 1

        assert self._degree >= 2 and self._degree <= 5, 'Degree of polynom must be between 3 and 5!'

        self._roots = []
        mpmath.mp.dps = dps

        #flip array
        coefficients_flipped = [self._coefficients[self._degree  - i] for i in range(0, self._degree + 1)]
        roots = polyroots(coefficients_flipped, extraprec=40, maxsteps=500, error = True)
        
        for root in roots[0]:
            self._roots.append(root)

        if self._degree == 5:
            self._permutation_group = Sym5()
        elif self._degree == 4:
            self._permutation_group = Sym4()
        elif self._degree == 3:
            self._permutation_group = Sym3()

        self._permutation_group = SymmetricGroup(self._degree)
        self._galois_resolvent = None
        self._sym_galois = None

    def __call__(self, x):
        return poly.polyval(x, self._coefficients)

    def galois_resolvent(self, max_iterations = MAX_ITERATIONS):
        iteration = 0
        found_galois_resolvente = False

        while not found_galois_resolvente and iteration < max_iterations:
            iteration += 1
            
            m  = [randint(-2, +2) for i in range(0, len(self._roots))]
            galois_res = GaloisResolvent(m, self._roots)
            number_equal = 0

            for i, per_i in enumerate(self._permutation_group._elements):
                galois_res_i = galois_res.permutate(per_i)

                for j, per_j in enumerate(self._permutation_group._elements):
                    if j > i:
                        galois_res_j = galois_res.permutate(per_j)

                        if galois_res_i.Value == galois_res_j.Value:
                            number_equal += 1
            
            if number_equal == 0:
                found_galois_resolvente = True
                self._galois_resolvent = galois_res

        if found_galois_resolvente == True:
            return 'Found Galois resolvent after {} iterations.'.format(str(iteration)), galois_res
        else:
            return 'Maximal number of iterations reached ({})! Could not find Galois resolvent.'.format(str(max_iterations))

    def sym_galois_poly(self):
        if self._galois_resolvent == None:
            _, self._galois_resolvent = self.galois_resolvent()
        
        sym_galois_pol = P([Decimal(1)])
        for i, per_i in enumerate(self._permutation_group._elements):
            galois_res = self._galois_resolvent.permutate(per_i)
            sym_galois_pol *= P([-galois_res.Value, Decimal(1)])

        self._sym_galois_pol = P([int(c) for c in sym_galois_pol.coef])
        return self._sym_galois_pol

    def get_factors(self, sym):
        subgroups = sym.subgroups()
        integer_polynoms = {}

        if self._galois_resolvent == None:
            _, self._galois_resolvent = self.galois_resolvent()

        for key, per_groups in subgroups.items():
            for per_group in per_groups:
                p = P([1.0])
                for per in per_group._elements:
                    galois_res = self._galois_resolvent.permutate(per).Value
                    p *= P([galois_res, 1.0])

                has_integer_coefficients, p_int = poly_over_integers(p)

                if has_integer_coefficients and is_factor(self._polynom, p_int):
                    if key not in integer_polynoms:
                        integer_polynoms[key] = [p_int, per_group._elements]
                    else:
                        integer_polynoms[key].append([p_int, per_group._elements])
        
        return integer_polynoms