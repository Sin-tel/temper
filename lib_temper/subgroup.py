# subgroups of the rational numbers

from fractions import Fraction
import numpy as np

from .util_types import Subgroup, IntMat
from .primes import primes
from .interval import *


# prime limit
def p_limit(limit: int) -> Subgroup:
    assert limit > 1
    assert limit < 7920

    return [i for i in primes if i <= limit]


# get basis for rational subgroup in term of prime expansion
def get_subgroup_basis(subgroup: Subgroup) -> tuple[IntMat, Subgroup]:
    expanded = prime_expansion(subgroup)
    s_basis = []
    for k in subgroup:
        s_basis.append(factors(k, expanded))
    s_basis = np.hstack(s_basis)

    # s_basis_new = temper.hnf(s_basis.T, remove_zeros=True).T

    # if redundant, normalize
    # if s_basis_new.shape[1] < s_basis.shape[1]:
    #   s_basis = s_basis_new

    # # if diagonal, normalize
    # if np.all(s_basis_new == np.diag(np.diagonal(s_basis_new))):
    #   s_basis = s_basis_new

    return s_basis, expanded


# expansion for fractional subgroup
def prime_expansion(subgroup: Subgroup):
    s: set[int] = set()

    for p in subgroup:
        s = s.union(prime_factors(p))

    return sorted(list(s))


# get the set of prime factors for a fraction. doesnt give multiplicity, just the primes
def prime_factors(frac: Fraction) -> int:
    s: set[int] = set()

    n, d = frac.as_integer_ratio()

    for p in primes:
        while n % p == 0:
            n //= p
            s.add(p)
        while d % p == 0:
            d //= p
            s.add(p)
        if n == 1 and d == 1:
            break

    assert n == 1, "Prime decomposition failed."
    assert d == 1, "Prime decomposition failed."

    return s


# get the rational subgroup from the prime expansion + basis
def get_subgroup(s_basis: IntMat, expanded: Subgroup) -> Subgroup:
    assert type(expanded) is list, "Subgroup must be a list."
    assert all(isinstance(x, int) for x in expanded)

    s = []
    for b in s_basis.T:
        s.append(ratio(b, expanded))
    return s
