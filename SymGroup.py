import numpy as np
import json
import os
from os import path

from sympy.combinatorics.named_groups import SymmetricGroup
from sympy.combinatorics import PermutationGroup
from sympy.combinatorics.permutations import Permutation
from helpers.math_helpers import *


class SymGroup:

    def __init__(self, n):
        self._sym_n = SymmetricGroup(n)
        self._elements = self._sym_n._elements
        self._permutation_group = PermutationGroup(self._sym_n._elements)
        self._path_serialization = os.path.join(os.path.dirname(__file__), 'symmetric_groups/sym_{0}.json'.format(n))
        self._subgroups = self.subgroups()
        self._trivial_group = PermutationGroup(Permutation([i for i in range(0, n)]))    

        with open(self._path_serialization, 'w') as write_file:
            json.dump(self._subgroups, write_file)

    @property
    def elements(self):
        return self._elements

    def find_maximal_normal_subgroup(self, subgroup):
        normal_subgroups = self.list_normal_subgroups(subgroup)
        max_order = 0
        max_normal_subgroup = self._trivial_group

        for normal_subgroup in normal_subgroups:
            order = len(normal_subgroup._elements)
            if order > max_order and order < len(subgroup._elements):
                max_order = order
                max_normal_subgroup = normal_subgroup

        return max_normal_subgroup
        
    def list_normal_subgroups(self, subgroup):        
        normal_subgroups = []
        subgroups_of = self.subgroups_of(subgroup)
        
        for subgroup_of in subgroups_of:
            normal_closure = subgroup.normal_closure(subgroup_of)

            if normal_closure == subgroup_of and not subgroup_of in normal_subgroups:
                normal_subgroups.append(subgroup_of)
                
        return normal_subgroups

    def subgroups_of(self, subgroup):
        all_subgroups = self.subgroups()
        subgroups_of = []
        
        for key, value in all_subgroups.items():
            for group in value:
                if group.is_subgroup(subgroup) and not group in subgroups_of:
                    subgroups_of.append(group)
                
        return subgroups_of

    def is_solvable_algebraically(self, subgroup):
        descending_groups = [subgroup]
        indices = []

        if subgroup.commutator(subgroup, ubgroup) == subgroup:
            return False, descending_groups, indices

        normal_group = subgroup
        normal_subgroup = self.find_maximal_normal_subgroup(normal_group)
        index = len(normal_group._elements) / len(normal_subgroup._elements)
        prime = is_prime(index)
        
        while prime:
            descending_groups.append(normal_subgroup)
            indices.append(index)

            if normal_subgroup.commutator(normal_subgroug, normal_subgroup) == normal_subgroup and not normal_subgroup == self._trivial_group:
                return False, descending_groups, indices

            normal_group = normal_subgroup
            normal_subgroup = self.find_maximal_normal_subgroup(normal_subgroup)
            index = len(normal_group._elements) / len(normal_subgroup._elements)
            prime = is_prime(index)

        if normal_subgroup == self._trivial_group:
            return True, descending_groups, indices
        else:
            return False, descending_groups, indices

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
                        auto_class_name = auto_class_name.replace('group_', '')
                        conj_subgroups = conjugated_subgroups(group_repr, auto_class_name, self._sym_n)

                        subgroups[auto_class_name] = {  'number conjugated subgroups' : group_repr.degree,
                                                        'group order' : len(group_repr._elements),    
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
        if not self._subgroups == None:
            return self._subgroups

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

