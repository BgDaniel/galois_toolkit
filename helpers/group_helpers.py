from sympy.combinatorics.named_groups import SymmetricGroup
from sympy.combinatorics import PermutationGroup
from sympy.combinatorics.permutations import Permutation
from helpers.math_helpers import *
class GroupTower:

    @property
    def DescendingGroups(self):
        return self._descending_groups

    def __init__(self, descending_groups, indices):
        self._descending_groups = descending_groups
        self._indices = indices

    def _iter__(self):
        pass
        
def is_solvable(group):
    descending_groups = [group]
    commutator = group.commutator(group, group)
    indices = []

    if commutator == group:
        return False, GroupTower(descending_groups)

    while len(commutator._elements) != 1:
        descending_groups.append(commutator)
        commutator_next = group.commutator(commutator, commutator)

        if commutator_next == commutator:
            return False, GroupTower(descending_groups)
        else:
            commutator = commutator_next

    # check whether G_i modulo G_i-1 is cyclic of prime order
    for i, group in enumerate(descending_groups):
        if i == len(descending_groups) - 1:
            continue
        index = len(group._elements) / len(descending_groups[i + 1]._elements)
        if not is_prime(index):
            return False, GroupTower(descending_groups, indices)
        else:
            indices.append(index)
    
    return True, GroupTower(descending_groups, indices)

def inv(permutation):
    array_inv = numpy.empty(len(permutation.array_form))

    for i, k in enumerate(permutation.array_form):
        array_inv[k] = i
    
    return Permutation(array_inv)

def mult(permutation_a, permutation_b):
    return Permutation([permutation_a[j] for i, j in enumerate(permutation_b.array_forms)])

class QuotientGroup:
    def __init__(self, group, sub_group):
        self._group = group
        self._sub_group = sub_group

    def equivalent(self, a, b):
        return mult(a, inv(b)) in self._sub_group

    def list_representatives(self):
        representatives = []

        for element in self._group._elements:
            for repr in representatives:
                if not equivalent(element, repr):
                    representatives.append(element)

        return representatives


    



  


