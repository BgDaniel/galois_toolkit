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

def power_sets(n, k):
    power_set = []
    assert k <= n, 'Cardinality of subset must be smaller than cardinality of super set!'
    if k == 0:
        return [[0 for j in range(0, n)]]
    elif k == n:
        return [[1 for j in range(0, n)]]
    elif k == 1:
        for i in range(0, n):
            subset = [1 if i == j else 0 for j in range(0, n)]
            power_set.append(subset)
        return power_set
    elif k == n - 1:
        for i in range(0, n):
            subset = [0 if i == j else 1 for j in range(0, n)]
            power_set.append(subset)
        return power_set
    else:
        for s_0 in power_sets(n - 1, k - 1):
            s_0_ext = s_0.copy()
            s_0_ext.append(1)
            power_set.append(s_0_ext)

        for s_1 in power_sets(n - 1, k):
            s_1_ext = s_1.copy()
            s_1_ext.append(0)
            power_set.append(s_1_ext)

        return power_set

def all_power_sets(n):
    all_power_sets = []

    for k in range(0, n + 1):
        power_sets_k = power_sets(n, k)
        for power_set_k in power_sets_k:
            all_power_sets.append(power_set_k)

    assert len(all_power_sets) == 2 ** n, 'Something went wrong in all_power_sets()!'
    return all_power_sets


