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

def poly_over_integers(p, precision = 10e-15):
    assert type(p) == PolyEquation or type(p) == P, 'Wrong input type for poly_over_integers()!'
    
    coefficients = []    
    
    if type(p) == PolyEquation:
        coefficients = p.Coefficients
    elif type(p) == P:
        coefficients = p.coef

    c_rounded = []
    for i, c in enumerate(coefficients):
        if abs(round(c) - c) > precision:
            return False, None
        else:
            c_rounded.append(int(round(c)))

    return True, P(c_rounded) 

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

        self._galois_resolvent = None
        self._sym_galois = None

    def from_polynom(polynom):
        return PolyEquation(polynom.coef)

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

            for i, per_i in enumerate(self._permutation_group.elements):
                galois_res_i = galois_res.permutate(per_i)

                for j, per_j in enumerate(self._permutation_group.elements):
                    if j > i:
                        galois_res_j = galois_res.permutate(per_j)

                        if galois_res_i.Value == galois_res_j.Value:
                            number_equal += 1
            
            if number_equal == 0:
                found_galois_resolvente = True
                self._galois_resolvent = galois_res

        if found_galois_resolvente == True:
            self.sym_galois_poly()
            return 'Found Galois resolvent after {} iterations.'.format(str(iteration)), galois_res
        else:
            return 'Maximal number of iterations reached ({})! Could not find Galois resolvent.'.format(str(max_iterations))

    def sym_galois_poly(self):
        if self._galois_resolvent == None:
            _, self._galois_resolvent = self.galois_resolvent()
        
        sym_galois_pol = P([1.0])
        for i, per_i in enumerate(self._permutation_group._elements):
            galois_res = self._galois_resolvent.permutate(per_i)
            sym_galois_pol *= P([-galois_res.Value, 1.0])

        self._sym_galois_pol = P([int(c) for c in sym_galois_pol.coef])
        return self._sym_galois_pol

    def determine_galois_group(self):
        integer_polynoms = {}
        galois_group = None

        if self._galois_resolvent == None:
            _, self._galois_resolvent = self.galois_resolvent()

        auto_classes = self._permutation_group.automorphism_classes()

        for key, value in auto_classes.items():
            p = P([1.0])
            
            for permutation in value._elements:
                galois_res = self._galois_resolvent.permutate(permutation).Value
                p *= P([galois_res, 1.0])

            p_over_int, p_int = poly_over_integers(p)

            if p_over_int:
                _, remainder = divmod(self._sym_galois_pol, p_int)

                if poly_over_integers(remainder):
                    integer_polynoms[key] = {   'polynom'       : str(p_int), 
                                                'group'         : key,
                                                'group order'   : len(value._elements),
                                                'elements'      : [per.array_form for per in value._elements],                                                 
                                            }

                    if galois_group == None:
                        galois_group = value
                    else:
                        if len(value._elements) < len(galois_group._elements):
                            galois_group = value

        return galois_group, integer_polynoms

        




