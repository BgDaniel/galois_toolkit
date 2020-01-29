import os
import json
from PolyEquation import *
from algebra.fields.algebraic_fields import *
from helpers.group_helpers import *
from symmetric_groups.Sym5 import *
from SymGroup import SymGroup



sym5 = Sym5()
wholeSym5 = sym5.group_Sym5()
descending_groups, indices = sym5.is_solvable_algebraically(wholeSym5)


coefficients = [1, 0, -10, 0, 1]
polynom = PolyEquation(coefficients, 75, 5000, -4, +4)

algebraic_field = AlgebraicField(polynom)
is_solvable, galois_group, group_tower = algebraic_field.determine_galois_group()


groups = group_tower.DescendingGroups

a = groups[0]
b = groups[1]

quotientGroup = QuotientGroup(b, a)



print(quotientGroup.list_representatives())










