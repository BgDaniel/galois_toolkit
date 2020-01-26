def ext_tuples(tuples, m_min, m_max):
    tuples_ext = []

    if len(tuples) == 0:
        for i in range(m_min, m_max):
            tuples_ext.append([i])

        return tuples_ext

    for tple in tuples:        
        for i in range(m_min, m_max):
            tple_cpy = list(tple)
            tple_cpy.append(i)
            tuples_ext.append(tple_cpy)

    return tuples_ext

def tuples(dim, m_min, m_max):
    tuples = []

    for i in range(0, dim):
        tuples = ext_tuples(tuples, m_min, m_max)

    return tuples

class Set:
    def __init__(self, order):
        self._order = order

    def power_set(self, k):
        assert k <= self._order, 'Wrong input!'
        power_set = []
        if self._order == 0:
            return [[0 for i in range(0, self._order)]]
        elif k == self._order:
            return [[1 for i in range(0, self._order)]]
        elif k == 1:
            for i in range(0, self._order):
                subset = [1 if i == j else 0 for j in range(0, self._order)]
                power_set.append(subset)
            return power_set
        else:
            power_set_inf = Set(self._order - 1)
            for s_0 in power_set_inf.power_set(k - 1):
                s_0_ext = s_0.copy()
                s_0_ext.append(1)
                power_set.append(s_0_ext)

            for s_1 in power_set_inf.power_set(k):
                s_1_ext = s_1.copy()
                s_1_ext.append(0)
                power_set.append(s_1_ext)

            return power_set