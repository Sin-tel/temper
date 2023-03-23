## interval arithmetic

import numpy as np
from fractions import Fraction


# prime vector to ratio. subgroup may be rational
def ratio(v, subgroup):
    v = v.flatten().tolist()
    p = 1
    q = 1
    for i, n in enumerate(v):
        if n > 0:
            p = p * subgroup[i] ** n
        elif n < 0:
            q = q * subgroup[i] ** (-n)

    return Fraction(p, q)


# fraction to prime vector (aka prime factorization)
# only works for integer subgroups (doesn't check for (co)primality in the basis)
def factors(fr, subgroup):
    assert type(subgroup) is list, "Subgroup must be a list."
    assert all(isinstance(x, int) for x in subgroup)

    factors = np.zeros((len(subgroup), 1), dtype=np.int64)

    if type(fr) is tuple:
        p = fr[0]
        q = fr[1]
    else:
        frac = Fraction(fr)
        p = frac.numerator
        q = frac.denominator

    assert p > 0
    assert q > 0

    for i, f in enumerate(subgroup):
        while p % f == 0:
            p //= f
            factors[i, 0] += 1
        while q % f == 0:
            q //= f
            factors[i, 0] -= 1

    # print(p,q, flush = True)
    assert p == 1 and q == 1, "Decomposition not in subgroup."

    return factors


def factors_unchecked(fr, subgroup):
    assert type(subgroup) is list, "Subgroup must be a list."
    assert all(isinstance(x, int) for x in subgroup)

    factors = np.zeros((len(subgroup), 1), dtype=np.int64)

    if type(fr) is tuple:
        p = fr[0]
        q = fr[1]
    else:
        frac = Fraction(fr)
        p = frac.numerator
        q = frac.denominator

    assert p > 0
    assert q > 0

    for i, f in enumerate(subgroup):
        while p % f == 0:
            p //= f
            factors[i, 0] += 1
        while q % f == 0:
            q //= f
            factors[i, 0] -= 1

    # print(p,q, flush = True)
    # assert p == 1 and q == 1, "Decomposition not in subgroup."
    if p != 1 or q != 1:
        return None

    return factors


# make a prime vector positive
def make_positive(v, subgroup):
    if log_interval(v, subgroup) < 0:
        return -v
    else:
        return v


# prime vector span in octaves
def log_interval(v, subgroup):
    vn = v.flatten()
    logs = log_subgroup(subgroup)

    return np.sum(vn * logs)


# subgroup in octaves as numpy array
def log_subgroup(subgroup):
    assert type(subgroup) is list, "Subgroup must be a list."

    return np.log2(np.array(subgroup).astype(np.double))
