from PolyEquation import *


coefficients = [-2, 8, -4, -4, 1]
polynom = PolyEquation(coefficients, 50, 5000, -2, +2)
galois_resolvent = polynom.galois_resolvent()
print(polynom.determine_galois_group())

coefficients = [-1, -2, 1, 1]
polynom = PolyEquation(coefficients, 50, 5000, -2, +2)
galois_resolvent = polynom.galois_resolvent()
print(polynom.determine_galois_group())

coefficients = [-1, -1, 0, 0, 0, 1]
polynom = PolyEquation(coefficients, 100, 5000, -5, +5)
galois_resolvent = polynom.galois_resolvent()
print(polynom.determine_galois_group())

coefficients = [1, 3, -3, -4, 1, 1]
polynom = PolyEquation(coefficients, 100, 5000, -5, +5)
galois_resolvent = polynom.galois_resolvent()
print(polynom.determine_galois_group())

coefficients = [1, 3, -3, -4, 1, 1]
polynom = PolyEquation(coefficients, 100, 5000, -5, +5)
galois_resolvent = polynom.galois_resolvent()
print(polynom.determine_galois_group())

