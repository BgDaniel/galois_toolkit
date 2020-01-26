import os
import json
from PolyEquation import *

list_coefficients = tuples(3, -2, 2)
result_file = os.path.join(os.path.dirname(__file__), 'results/results.json')

results = []

with open(result_file, 'w') as write_file:
    json.dump(results, write_file)

for coefficients in list_coefficients:
    coefficients.append(1)
    polynom = PolyEquation(coefficients, 100, 5000, -7, +7)
    group_class, solvable, galois_group, _ = polynom.determine_galois_group()

    if group_class == None:
        continue

    result = { 'coefficients'   : coefficients,
               'group class'    : group_class,
               'solvable'       : solvable,
               'galois group'   : galois_group 
             }
    
    with open(result_file, 'r') as read_file:
        results_str = read_file.read()
        results = json.loads(results_str)

    results.append(result)

    with open(result_file, 'w') as write_file:
        json.dump(results, write_file)




