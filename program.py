from PolyEquation import *


coefficients = [-2, 8, -4, -4, 1]
polynom = PolyEquation(coefficients)
galois_resolvent = polynom.galois_resolvent()
polynom.determine_galois_group()

