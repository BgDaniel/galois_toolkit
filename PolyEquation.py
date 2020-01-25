import numpy.polynomial.polynomial as poly
from numpy.polynomial.polynomial import Polynomial as P
from sympy.combinatorics.named_groups import SymmetricGroup
import mpmath
from mpmath import *
from random import randint

from GaloisResolvent import *
from symmetric_groups.Sym5 import Sym5
from symmetric_groups.Sym4 import Sym4
from symmetric_groups.Sym3 import Sym3
from helpers.math_helpers import *
from helpers.group_helpers import *


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
        if type(c) == mpc:
            if abs(nint(c.imag) - .0) > precision:
                return False, None
            if abs(nint(c.real) - c.real) > precision:
                return False, None
            else:
                c_rounded.append(nint(c.real))
        else:
            if abs(nint(c) - c) > precision:
                return False, None
            else:
                c_rounded.append(nint(c))

    return True, P(c_rounded) 

def poly_to_str(p):
    assert type(p) == PolyEquation or type(p) == P, 'Wrong input type for poly_over_integers()!'
    p_str = ""

    coefficients = []    
    
    if type(p) == PolyEquation:
        coefficients = p.Coefficients
    elif type(p) == P:
        coefficients = p.coef

    for i, c in enumerate(coefficients):
        if c == 0:
            continue
        if i == 0:
            if c < 0:
                p_str = "- {0} ".format(-c) + p_str
            else:
                p_str = "+ {0} ".format(c) + p_str
        elif i == 1:
            if c < 0:
                p_str = "- {0} * X ".format(-c) + p_str
            else:
                p_str = "+ {0} * X ".format(c) + p_str
        else:
            if c < 0:
                if c == 1.0:
                    p_str = "- X^{1} ".format(-c, i) + p_str
                else:
                    p_str = "- {0} * X^{1} ".format(-c, i) + p_str
            else:
                if c == 1.0:
                    if i == len(coefficients) - 1:
                        p_str = "X^{1} ".format(-c, i) + p_str
                    else:
                        p_str = "+ X^{1} ".format(-c, i) + p_str
                else:
                    p_str = "+ {0} * X^{1} ".format(c, i) + p_str

    return p_str

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

    def __init__(self, coefficients, dps, max_iterations, m_min, m_max):
        self._coefficients = coefficients
        self._polynom = P(coefficients) 
        self._degree = len(coefficients) - 1
        self._max_iterations = max_iterations
        self._m_min = m_min
        self._m_max = m_max

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
        self._integer_polynoms = None

    def from_polynom(polynom):
        return PolyEquation(polynom.coef)

    def __call__(self, x):
        return poly.polyval(x, self._coefficients)

    def galois_resolvent(self):
        iteration = 0
        number_equal = - 1

        for m in tuples(self._degree, self._m_min, self._m_max):
            iteration += 1
            
            galois_res = GaloisResolvent(m, self._roots)
            number_equal = 0

            for i, per_i in enumerate(self._permutation_group.elements):
                if number_equal > 0:
                    continue
                
                galois_res_i = galois_res.permutate(per_i)

                for j, per_j in enumerate(self._permutation_group.elements):
                    if number_equal > 0:
                        continue

                    if j > i:
                        galois_res_j = galois_res.permutate(per_j)

                        if galois_res_i.Value == galois_res_j.Value:
                            number_equal += 1
            
            if number_equal == 0:
                self._galois_resolvent = galois_res
                self.sym_galois_poly()
                return 'Found Galois resolvent after {} iterations.'.format(str(iteration)), galois_res
            elif self._max_iterations == iteration:
                return 'Maximal number of iterations reached ({})! Could not find Galois resolvent.'.format(str(max_iterations))

    def sym_galois_poly(self):
        if self._galois_resolvent == None:
            _, self._galois_resolvent = self.galois_resolvent()
        
        sym_galois_pol = P([1.0])
        for i, per_i in enumerate(self._permutation_group._elements):
            galois_res = self._galois_resolvent.permutate(per_i)
            sym_galois_pol *= P([-galois_res.Value, 1.0])

        over_int, sym_galois_pol_int = poly_over_integers(sym_galois_pol)

        if over_int:
            self._sym_galois_pol = P([c for c in sym_galois_pol.coef])
            return self._sym_galois_pol
        else:
            raise Exception('Nummerically calculated symmetric Galois polynom has no integer coefficients!')

    def determine_galois_group(self):
        self._integer_polynoms = {}
        galois_group = None
        galois_polynom = None
        galois_name = None
        galois_is_solvable = False

        if self._galois_resolvent == None:
            _, self._galois_resolvent = self.galois_resolvent()

        subgroups = self._permutation_group.subgroups()

        for key, value in subgroups.items():
            for subgroup in value:
                p = P([1.0])
            
                for permutation in subgroup._elements:
                    galois_res = self._galois_resolvent.permutate(permutation).Value
                    p *= P([galois_res, 1.0])

                p_over_int, p_int = poly_over_integers(p)

                if p_over_int:
                    _, remainder = divmod(self._sym_galois_pol, p_int)

                    if poly_over_integers(remainder):
                        is_solv, _ = is_solvable(subgroup)
                        self._integer_polynoms[key] =   {   
                                                        'polynom'       : str(p_int), 
                                                        'group'         : key,
                                                        'group order'   : len(subgroup._elements),
                                                        'solvable'      : is_solv,
                                                        'elements'      : [per.array_form for per in subgroup._elements],                                                 
                                                    }

                        if galois_group == None:
                            galois_group = subgroup
                            galois_polynom = p_int
                            galois_name = key
                            galois_is_solvable = is_solv
                        else:
                            if len(subgroup._elements) < len(galois_group._elements):
                                galois_group = subgroup
                                galois_polynom = p_int
                                galois_name = key
                                galois_is_solvable = is_solv

        return galois_name, "solvable: {0}".format(galois_is_solvable), [per.array_form for per in galois_group._elements], poly_to_str(galois_polynom)

        




