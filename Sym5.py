from sympy.combinatorics import PermutationGroup
from SymGroup import SymGroup
from sympy.combinatorics.permutations import Permutation

class Sym5(SymGroup):
    def __init__(self):
        SymGroup.__init__(self, 5)

    '''https://groupprops.subwiki.org/wiki/Symmetric_group:S5'''

    '''trivial subgroup'''
    def group_Trivial(self):
        return PermutationGroup([Permutation([0, 1, 2, 3])])

    def group_S2(self):
        return PermutationGroup([Permutation([1, 0, 2, 3])])

    def group_Dbl_Trans(self):
        return PermutationGroup([Permutation([1, 0, 3, 2])])

    def group_Z4(self):
        return PermutationGroup([Permutation([3, 0, 1, 2])])

    def group_Klein4_N(self):
        return PermutationGroup([Permutation([1, 0, 3, 2]), Permutation([2, 3, 0, 1]), Permutation([3, 2, 1, 0])])

    def group_Klein4(self):
        return PermutationGroup([Permutation([1, 0, 2, 3]), Permutation([0, 1, 3, 2])])

    def group_D8(self):
        return PermutationGroup([Permutation([3, 0, 1, 2]), Permutation([2, 1, 0, 3])])

    def group_A3(self):
        return PermutationGroup([Permutation([2, 0, 1, 3]), Permutation([1, 2, 0, 3])])

    def group_S3(self):
        return PermutationGroup([Permutation([2, 0, 1, 3]), Permutation([1, 0, 2, 3])])

    def group_A4(self):
        return PermutationGroup([Permutation([2, 0, 1, 3]), Permutation([1, 0, 3, 2])])




sym5 = Sym5()
auto_classes = subgroups = sym5.list_automorphism_class_of_subgroups()
print(auto_classes)