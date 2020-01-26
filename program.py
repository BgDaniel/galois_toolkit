import os
import json
from PolyEquation import *

coefficients = [4, 11, 2, -4, -2, 1]
polynom = PolyEquation(coefficients, 50, 5000, -2, +2)
print(polynom.is_irreducible())








