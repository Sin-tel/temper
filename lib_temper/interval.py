# interval arithmetic

from fractions import Fraction
import numpy as np
from .util_types import IntArray, IntVec, FloatVec, Subgroup, SubgroupInt, IntMat


# prime vector to ratio. subgroup may be rational
def ratio(v: IntVec, subgroup: Subgroup) -> Fraction:
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
def factors(fr: int | Fraction | tuple[int, int], subgroup: SubgroupInt) -> IntVec:
    f_vector = np.zeros((len(subgroup), 1), dtype=np.int64)

    if isinstance(fr, tuple):
        fr = Fraction(fr[0], fr[1])

    frac = Fraction(fr)
    p = frac.numerator
    q = frac.denominator

    assert p > 0
    assert q > 0

    for i, f in enumerate(subgroup):
        while p % f == 0:
            p //= f
            f_vector[i, 0] += 1
        while q % f == 0:
            q //= f
            f_vector[i, 0] -= 1

    assert p == 1 and q == 1, f"Decomposition of {fr} not in subgroup."

    return f_vector


# make a prime vector positive
def make_positive(v: IntVec, subgroup: Subgroup) -> IntVec:
    if log_interval(v, subgroup) < 0:
        return -v
    return v


# prime vector span in octaves
def log_interval(v: IntVec, subgroup: Subgroup) -> float:
    vn = v.flatten()
    logs = log_subgroup(subgroup)

    return np.sum(vn * logs)


# subgroup in octaves as numpy array
def log_subgroup(subgroup: Subgroup) -> FloatVec:
    return np.log2(np.array(subgroup).astype(np.float64))
