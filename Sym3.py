from sympy.combinatorics import PermutationGroup
from SymGroup import SymGroup
from sympy.combinatorics.permutations import Permutation

class Sym3(SymGroup):
    def __init__(self):
        SymGroup.__init__(self, 3)

    '''https://groupprops.subwiki.org/wiki/Subgroup_structure_of_symmetric_group:S3'''

    '''whole subgroup'''
    def group_Sym3(self):
        return self._sym_n

    '''trivial subgroup'''
    def group_Trivial(self):
        return PermutationGroup([Permutation([0, 1, 2])])

    '''S2 in S3'''
    def group_S2(self):
        return PermutationGroup([Permutation([1, 0, 2]), Permutation([0, 2, 1]), Permutation([2, 1, 0])])

    '''A3 in S3'''
    def group_A3(self):
        return PermutationGroup([Permutation([2, 0, 1]), Permutation([1, 2, 0])])
