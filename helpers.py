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