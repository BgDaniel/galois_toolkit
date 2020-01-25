import numpy as np
import numpy.polynomial.polynomial as poly
from numpy.polynomial.polynomial import Polynomial as P
from sympy.combinatorics.named_groups import SymmetricGroup


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


      




