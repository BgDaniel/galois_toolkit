import unittest
from PolyEquation import *

TEST_CASES = [
    { 'coefficients' : [-2, 8, -4, -4, 1], 'group' : 'D8', 'solvable' : 'solvable: True'  },
    { 'coefficients' : [-1, -2, 1, 1], 'group' : 'A3', 'solvable' : 'solvable: True'  },
    { 'coefficients' : [-1, -1, 0, 0, 0, 1], 'group' : 'Sym5', 'solvable' : 'solvable: False'  },
    { 'coefficients' : [1, 3, -3, -4, 1, 1], 'group' : 'Z5', 'solvable' : 'solvable: True'  },
    { 'coefficients' : [12, -5, 0, 0, 0, 1], 'group' : 'D10', 'solvable' : 'solvable: True'  }
]

class TestGaloisGroup(unittest.TestCase):
    def test_galois_group(self):
        for test_case in TEST_CASES:
            polynom = PolyEquation(test_case['coefficients'], 100, 5000, -5, +5)
            group_class, solvable, _, _ = polynom.determine_galois_group()
            self.assertEqual(group_class, test_case['group'])
            self.assertEqual(solvable, test_case['solvable'])
            
if __name__ == '__main__':
    unittest.main()