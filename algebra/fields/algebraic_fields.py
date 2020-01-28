from PolyEquation import PolyEquation

class AlgebraicField:
    @property
    def GaloisGroup(self):
        return self._galois_group

    @property
    def ISolvable(self):
        return self._is_solvable

    @property
    def GroupTower(self):
        return self._group_tower

    def __init__(self, minimal_polynomial):
        self._minimal_polynomial = minimal_polynomial
        self._roots = minimal_polynomial.Roots
        self._galois_group = None
        self._group_tower = None
        self._field_tower = None

    def determine_galois_group(self):
        self._is_solvable, self._galois_group, self._group_tower = self._minimal_polynomial.determine_galois_group()
        return self._is_solvable, self._galois_group, self._group_tower




