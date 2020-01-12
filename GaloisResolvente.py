import numpy as np
import numpy.polynomial.polynomial as poly
from numpy.polynomial.polynomial import Polynomial as P
from numpy.random import randint
from sympy.combinatorics.named_groups import SymmetricGroup
from sympy.combinatorics import PermutationGroup
import math
import decimal
from decimal import *
import mpmath
from mpmath import *
import mpmath as mp

from Sym4 import Sym4

MAX_ITERATIONS = 5000
PRECISION = 0.00000001



class GaloisResolvent:

    @property
    def Value(self):
        return self._value

    @property
    def Abs(self):
        return self._abs

    def __init__(self, m, roots):
        self._m = m
        self._roots = roots
        
        self._value = .0
        for i, r in enumerate(self._roots):
            self._value += r * self._m[i] 
        
        self._abs = abs(self._value)
        self._permutation_group = SymmetricGroup(len(roots))

    def permutate(self, permutation):
        per_array = permutation.array_form
        roots_permutated = [self._roots[per_array[i]] for i in range(0, len(self._roots))]
        return GaloisResolvent(np.array(self._m), roots_permutated)

    def orbit(self):
        orbit = []
        max_norm = .0

        for per in self._permutation_group._elements:
            galois_res_per = self.permutate(per)
            orbit.append(galois_res_per)

            abs = galois_res_per.Abs
            if abs > max_norm:
                max_norm = abs

        return max_norm, orbit

    def find(roots, max_norm = 2.0, max_iteration = 5000, m_min = - 250, m_max = + 250):
        found_solution = False
        iteration = 0

        while not found_solution and iteration < max_iteration:
            iteration += 1
            m  = randint(m_min, m_max, len(roots))
            galois_res = GaloisResolvent(m, roots)

            max_norm_orbit , orbit = galois_res.orbit()

            if max_norm_orbit < max_norm:
                found_solution = True

        if found_solution == True:
            print('Found Galois resolvent whose orbit is bounded by {0} after {1} iterations.'.format(str(max_norm), str(iteration)))
            return galois_res
        else:
            print('Maximal number of iterations reached ({0})! Could not find Galois resolvent whose orbit is bounded by {1}.'.format(str(max_iterations)), str(max_norm))





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

    def __init__(self, coefficients):
        self._coefficients = coefficients
        self._polynom = P(coefficients) 
        self._degree = len(coefficients) - 1
        self._roots = []
        mpmath.mp.dps = 35

        #flip array
        coefficients_flipped = [self._coefficients[self._degree  - i] for i in range(0, self._degree + 1)]
        roots = polyroots(coefficients_flipped, extraprec=40, maxsteps=500, error = True)
        
        for root in roots[0]:
            self._roots.append(root)

        self._degree = len(coefficients) - 1
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
            
            m  = randint(-2, +2, len(self._roots))
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

        

def poly_over_integers(p, precision = 10e-9):
    c_rounded = []
    for i, c in enumerate(p.coef):
        if abs(round(c) - c) > 10e-9:
            return False, None
        else:
            c_rounded.append(int(round(c)))

    return True, P(c_rounded) 

def is_factor(p, q):

    self._sym_galois
    s = divmod(p, q)
    return(s)


      

class SubgroupFinder:
    def __init__(self, permutation_group):
        self._permutation_group = permutation_group
        self._elements_super_group = self._permutation_group._elements
        self._order_super_group = len(self._permutation_group._elements)

    def list_subgroups(self, nb_generators):
        s = Set(self._order_super_group)
        power_set = s.power_set(nb_generators)

        subgroups = []
        for subset in power_set:
            generators = []
            
            for i, indicator in enumerate(subset):
                if indicator == 1:
                    generators.append(self._elements_super_group[i])
            
            permutation_subgroup = PermutationGroup(generators)

            if not permutation_subgroup in subgroups:
                subgroups.append(permutation_subgroup)

        return subgroups

    def list_all_subgroups(self):
        all_subgroups = []
        for nb_generators in range(1, self._order_super_group):
            for subgroup_of_order in self.list_subgroups(nb_generators):
                if not subgroup_of_order in all_subgroups:
                    all_subgroups.append(subgroup_of_order)

        sub_groups_struc = {}

        for subgroup in all_subgroups:
            order = subgroup._order
            key_order = 'order {}'.format(str(order))

            if not key_order in sub_groups_struc.keys():
                sub_groups_struc[key_order] = [subgroup]
            else:
                sub_groups_struc[key_order].append(subgroup)
            
        return all_subgroups



class Set:
    def __init__(self, order):
        self._order = order

    def power_set(self, k):
        assert k <= self._order, 'Wrong input!'
        power_set = []
        if self._order == 0:
            return [[0 for i in range(0, self._order)]]
        elif k == self._order:
            return [[1 for i in range(0, self._order)]]
        elif k == 1:
            for i in range(0, self._order):
                subset = [1 if i == j else 0 for j in range(0, self._order)]
                power_set.append(subset)
            return power_set
        else:
            power_set_inf = Set(self._order - 1)
            for s_0 in power_set_inf.power_set(k - 1):
                s_0_ext = s_0.copy()
                s_0_ext.append(1)
                power_set.append(s_0_ext)

            for s_1 in power_set_inf.power_set(k):
                s_1_ext = s_1.copy()
                s_1_ext.append(0)
                power_set.append(s_1_ext)

            return power_set





s = Set(5)
power_set = s.power_set(3)
print(power_set)


#sym_4 = SymmetricGroup(4)
#sub_group_finder = SubgroupFinder(sym_4)
#all_sub_groups = sub_group_finder.list_all_subgroups()





coefficients = [-2, 8, -4, -4, 1]
f = PolyEquation(coefficients)

sym4 = Sym4()
f.get_factors(sym4)

print(f.sym_galois_poly())

for root in f.Roots:
    print(f(root))
