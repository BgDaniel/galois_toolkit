import numpy as np
import json
import os
from os import path

from sympy.combinatorics.named_groups import SymmetricGroup
from sympy.combinatorics import PermutationGroup
from sympy.combinatorics.permutations import Permutation
from helpers.group_helpers import *

class SymGroup:

    def __init__(self, n):
        self._sym_n = SymmetricGroup(n)
        self._elements = self._sym_n._elements
        self._path_serialization = os.path.join(os.path.dirname(__file__), 'symmetric_groups/sym_{0}.json'.format(n))
        self._subgroups = self.list_subgroups()        

        with open(self._path_serialization, 'w') as write_file:
            json.dump(self._subgroups, write_file)

    @property
    def elements(self):
        return self._elements
        
    def list_subgroups(self):        
        if path.exists(self._path_serialization):
            with open(self._path_serialization, 'r') as read_file:
                subgroups_str = read_file.read()
                return json.loads(subgroups_str)
        else:
            subgroups = {}

            for auto_class_name in dir(self):
                if auto_class_name.startswith('group_'):
                    auto_class = getattr(self, auto_class_name)

                    if callable(auto_class):
                        group_repr = auto_class()
                        is_solv, group_tower = is_solvable(group_repr)
                        auto_class_name = auto_class_name.replace('group_', '')
                        conj_subgroups = conjugated_subgroups(group_repr, auto_class_name, self._sym_n)

                        subgroups[auto_class_name] = {  'number conjugated subgroups' : group_repr.degree,
                                                        'group order' : len(group_repr._elements),
                                                        'solvable': is_solv,     
                                                        'subgroups' : [[per.array_form for per in conj_subgroup._elements] for conj_subgroup in conj_subgroups]
                                                    }
                                        
            return subgroups

    def automorphism_classes(self):
        auto_classes = {}

        for auto_class_name in dir(self):
            if auto_class_name.startswith('group_'):
                auto_class = getattr(self, auto_class_name)

                if callable(auto_class):
                    group_repr = auto_class()
                    auto_class_name = auto_class_name.replace('group_', '')
                    auto_classes[auto_class_name] = group_repr

        return auto_classes

    def subgroups(self):
        subgroups = {}

        for auto_class_name in dir(self):
            if auto_class_name.startswith('group_'):
                auto_class = getattr(self, auto_class_name)

                if callable(auto_class):
                    group_repr = auto_class()
                    auto_class_name = auto_class_name.replace('group_', '')
                    subgroups[auto_class_name] = conjugated_subgroups(group_repr, auto_class_name, self._sym_n)

        return subgroups

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

    return conjugated_subgroups

