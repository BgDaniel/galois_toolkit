from sympy.combinatorics.named_groups import SymmetricGroup
from sympy.combinatorics import PermutationGroup


class SymGroup(PermutationGroup):
    def __init__(self, n):
        self._sym_n = SymmetricGroup(n)
        PermutationGroup.__init__(self._sym_n._elements)

    def list_subgroups(self):
        subgroups = { self._sym_n.order : [self._sym4] }

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

    def list_automorphism_class_of_subgroups(self):
        auto_classes = { 'whole group' : [self._sym_n] }

        for auto_class_name in dir(self):
            if auto_class_name.startswith('group_'):
                auto_class = getattr(self, auto_class_name)

                if callable(auto_class):
                    group_repr = auto_class()
                    auto_class_name = auto_class_name.replace('group_', '')
                    auto_classes[auto_class_name] = group_repr

        return auto_classes

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

