from sympy.combinatorics import PermutationGroup
from SymGroup import SymGroup
from sympy.combinatorics.permutations import Permutation

class Sym4(SymGroup):
    def __init__(self):
        SymGroup.__init__(self, 4)

    '''https://www.google.de/search?q=add+url+as+comment+python&ie=&oe='''

    '''trivial subgroup'''
    def group_Trivial(self):
        return PermutationGroup([Permutation([0, 1, 2, 3])])

    '''S2 in S4'''
    def group_S2(self):
        return PermutationGroup([Permutation([1, 0, 2, 3])])

    '''subgroup generated by double transposition in S4'''
    def group_Dbl_Trans(self):
        return PermutationGroup([Permutation([1, 0, 3, 2])])

    '''Z4 in S4'''
    def group_Z4(self):
        return PermutationGroup([Permutation([3, 0, 1, 2])])

    '''normal Klein four-subgroup of S4'''
    def group_Klein4_N(self):
        return PermutationGroup([Permutation([1, 0, 3, 2]), Permutation([2, 3, 0, 1]), Permutation([3, 2, 1, 0])])

    '''non-normal Klein four-subgroups of S4'''
    def group_Klein4(self):
        return PermutationGroup([Permutation([1, 0, 2, 3]), Permutation([0, 1, 3, 2])])

    '''D8 in S4'''
    def group_D8(self):
        return PermutationGroup([Permutation([3, 0, 1, 2]), Permutation([2, 1, 0, 3])])

    '''A3 in S4'''
    def group_A3(self):
        return PermutationGroup([Permutation([2, 0, 1, 3]), Permutation([1, 2, 0, 3])])

    '''S3 in S4'''
    def group_S3(self):
        return PermutationGroup([Permutation([2, 0, 1, 3]), Permutation([1, 0, 2, 3])])

    '''A4 in S4'''
    def group_A4(self):
        return PermutationGroup([Permutation([2, 0, 1, 3]), Permutation([1, 0, 3, 2])])


sym4 = Sym4()
auto_classes = subgroups = sym4.list_automorphism_class_of_subgroups()
print(auto_classes)