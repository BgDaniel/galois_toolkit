from sympy.combinatorics.named_groups import SymmetricGroup
from sympy.combinatorics import PermutationGroup
from sympy.combinatorics.permutations import Permutation

class GroupTower:

    def __init__(self, descending_groups):
        self._ascending_groups = descending_groups

    def _iter__(self):
        pass
        
def is_solvable(group):
    descending_groups = [group]
    commutator = group.commutator(group, group)

    if commutator == group:
        return False, GroupTower(descending_groups)

    while len(commutator._elements) != 1:
        descending_groups.append(commutator)
        commutator_next = commutator.commutator(commutator, commutator)

        if commutator_next == commutator:
            return False, GroupTower(descending_groups)
        else:
            commutator = commutator_next

    return True, GroupTower(descending_groups)


  


