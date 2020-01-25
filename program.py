from PolyEquation import *


coefficients = [-2, 8, -4, -4, 1]
f = PolyEquation(coefficients)
galois_resolvent = f.galois_resolvent()

sym4 = Sym4()
f.get_factors(sym4)

print(f.sym_galois_poly())

for root in f.Roots:
    print(f(root))