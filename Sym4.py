from sympy.combinatorics.named_groups import SymmetricGroup
from sympy.combinatorics import PermutationGroup
from enum import Enum
from sympy.combinatorics.permutations import Permutation
import numpy as np

def inv(p):
    ell = len(p.array_form)
    arr_inv = np.zeros(ell)
    
    for i, j in enumerate(p.array_form):
        arr_inv[j] = i

    return Permutation(arr_inv)

def conjugate_subgroup(g, subgroup):
    elements_conjugated = []
    
    for h in subgroup._elements:
        elements_conjugated.append(g * h * inv(g))

    return PermutationGroup(elements_conjugated)

def is_in(H, subgroups_list):
    for G in subgroups_list:
        if G == H:
            return True

    return False

def conjugated_subgroups(subgroup, name, permutation_group):
    conjugated_subgroups = [subgroup]

    for g in permutation_group._elements:
        subgroup_conjugated = conjugate_subgroup(g, subgroup)

        if not is_in(subgroup_conjugated, conjugated_subgroups):
            conjugated_subgroups.append(subgroup_conjugated)

        #if not subgroup_conjugated in conjugated_subgroups:
        #    conjugated_subgroups.append(subgroup_conjugated)

    print('The subgroup {} has {} conjugated subgroups.'.format(name, len(conjugated_subgroups)))
    return conjugated_subgroups

class Sym4(PermutationGroup):
    def __init__(self):
        self._sym4 = SymmetricGroup(4)
        PermutationGroup.__init__(self._sym4._elements)

    def subgroups(self):
        subgroups = { self._sym4.order : [self._sym4] }

        for auto_class_name in dir(self):
            if auto_class_name.startswith('group_'):
                auto_class = getattr(self, auto_class_name)

                if callable(auto_class):
                    group_repr = auto_class()
                    auto_class_name = auto_class_name.replace('group_', '')
                    conj_subgroups = conjugated_subgroups(group_repr, auto_class_name, self._sym4)
                    for conj_subgroup in conj_subgroups:
                        if not conj_subgroup in subgroups:
                            order = conj_subgroup.order()
                            if order in subgroups:
                                subgroups[order].append(conj_subgroup)
                            else:
                                subgroups[order] = [conj_subgroup]

        return subgroups

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


sym4 = Sym4()
subgroups = sym4.subgroups()