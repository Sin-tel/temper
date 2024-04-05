# collection of functions for dealing with regular temperaments

from typing import Optional, TypeVar

import numpy as np

from . import olll
from . import diophantine
from .hnf_bigint import hnf_bigint

from .util_types import IntMat, FloatMat
from .subgroup import *
from .optimize import *
from .combo import comboBySum


# Find the hermite normal form of M
def hnf(M: IntMat, remove_zeros: bool = False, transformation: bool = False) -> IntMat:
    assert M.ndim == 2

    try:
        if transformation:
            solution = diophantine.lllhermite(M.astype(np.int64))
            res = np.array(solution[0]).astype(np.int64)
            unimod = np.array(solution[1]).astype(np.int64)
        else:
            res = hnf_bigint(M)
    except OverflowError as e:
        raise OverflowError("Your numbers are too big!") from e

    if remove_zeros:
        idx = np.argwhere(np.all(res[:] == 0, axis=1))
        res = np.delete(res, idx, axis=0)

    if transformation:
        return res, unimod  # type: ignore[return-value]
    return res


# Find the left kernel (nullspace) of M.
# Adjoin an identity block matrix and solve for HNF.
# This is equivalent to the highschool maths method,
# but using HNF instead of Gaussian elimination.
def kernel(M: IntMat) -> IntMat:
    assert M.ndim == 2
    r, d = M.shape

    M = np.vstack([M, np.eye(d, dtype=np.int64)])
    K = hnf(M.T).T
    K = K[r::, r::]

    return K


# Find the right kernel (nullspace) of M
def cokernel(M: IntMat) -> IntMat:
    assert M.ndim == 2
    return kernel(M.T).T


T = TypeVar("T", bound=IntMat | FloatMat)


# LLL reduction
# Tries to find the smallest basis vectors in some lattice
# Operates on columns
# M: Map
# W: Weight matrix (metric)
def LLL(M: T, W: Optional[FloatMat] = None, delta: float = 0.99) -> T:
    if W is None:
        W = np.eye(M.shape[0])
    res = olll.reduction(np.copy(M).T, delta=delta, W=W).T

    # sort them by complexity
    # actually, this might be redundant.
    c_list = list(res.T)
    c_list.sort(key=lambda c: np.dot(c, W @ c))

    return np.array(c_list).T


# Flips M along diagonal
def antitranspose(M):
    return np.flipud(np.fliplr((M.T)))


# Finds the Hermite normal form and 'defactors' it.
# Defactoring is also known as saturation.
# This removes torsion from the map.
# Algorithm as described by:
#
# Clément Pernet and William Stein.
# Fast Computation of HNF of Random Integer Matrices.
# Journal of Number Theory.
# https://doi.org/10.1016/j.jnt.2010.01.017
# See section 8.
def defactored_hnf(M):
    r, d = M.shape

    S = np.linalg.inv(hnf(M.T)[:r].T)

    assert np.allclose(S @ M, np.round(S @ M))

    D = np.round(S @ M).astype(np.int64)

    return hnf(D)


# exact integer determinant using Bareiss algorithm
# modified slightly from:
## https://stackoverflow.com/questions/66192894/precise-determinant-of-integer-nxn-matrix
def integer_det(M):
    M = np.copy(M)  # make a copy to keep original M unmodified

    N, sign, prev = len(M), 1, 1
    for i in range(N - 1):
        if M[i, i] == 0:  # swap with another row having nonzero i's elem
            swapto = next((j for j in range(i + 1, N) if M[j, i] != 0), None)
            if swapto is None:
                return 0  # all M[*][i] are zero => zero determinant
            ## swap rows
            M[[i, swapto]] = M[[swapto, i]]
            sign *= -1
        for j in range(i + 1, N):
            for k in range(i + 1, N):
                assert (M[j, k] * M[i, i] - M[j, i] * M[i, k]) % prev == 0
                M[j, k] = (M[j, k] * M[i, i] - M[j, i] * M[i, k]) // prev
        prev = M[i, i]
    return sign * M[-1, -1]


# Order of factorization.
# For a saturated basis this is 1.
def factor_order(M):
    r = M.shape[0]
    return integer_det(hnf(M.T)[:r].T)


# Canonical maps
# This is just the defactored HNF,
# but for comma bases we do the antitranspose sandwich.
def canonical(M):
    assert M.ndim == 2
    r, d = M.shape
    if r > d:
        # comma basis
        return antitranspose(defactored_hnf(antitranspose(M)))
    # mapping
    return defactored_hnf(M)


# Solve AX = B in the integers
# for the method used, see https://github.com/tclose/Diophantine/blob/master/algorithm.pdf
def solve_diophantine(A, B):
    B = np.atleast_2d(B)
    assert A.shape[0] == B.shape[0]
    aug = np.block(
        [
            [A.T, np.zeros((A.shape[1], B.shape[1]), dtype=np.int64)],
            [B.T, np.eye(B.shape[1], dtype=np.int64)],
        ]
    )

    r, d = A.shape

    nullity = d - r

    # somehow we can solve the system even if the nullity is -1
    # aka the kernel is trivial
    # should double check when this actually works
    if nullity <= 0:
        nullity = 0

    _, U = hnf(aug, transformation=True)

    p2 = U.shape[0] - nullity
    p1 = p2 - B.shape[1]

    sol = -U[p1:p2, : A.shape[1]].T

    # Check that the solution actually works.
    # Probably the easiest way to guarantee this routine works correctly.
    assert np.all(A @ sol == B), "Could not solve system"

    return sol


def preimage(M):
    # I have no idea how I came up with this but it seems to work

    r, d = M.shape

    B = np.block([[M.T, np.eye(d, dtype=np.int64)]])
    H = hnf(B)
    sol = H[0:r, r:].T

    assert np.all(M @ sol == np.eye(r, dtype=np.int64)), "Could not solve system"

    return sol


# simplify list of intervals wrt a comma basis.
# The comma basis should be in reduced LLL form for this to work properly.
# TODO: figure out if the weight matrix here actually helps
def simplify(intervals, commas, W=None):
    # system tries to find vectors approximately equal to zero
    # approximate solution in floating point,
    # then round to get integer solution

    if W is None:
        W = np.eye(commas.shape[0])

    # c_inv = np.linalg.pinv(commas)
    c_inv = np.linalg.inv(commas.T @ W @ commas) @ commas.T
    s_vecs = np.round(c_inv @ W @ intervals).astype(np.int64)
    simpl = intervals - commas @ s_vecs
    return simpl


# Patent n-edo map
# Just crudely rounds all the log primes, multiplied by n
def patent_map(edo_num, subgroup):
    log_s = log_subgroup(subgroup)
    log_s = log_s / log_s[0]  # fix equave

    # floor(x+0.5) rounds more predictably (downwards on .5)
    M = np.floor(edo_num * log_s + 0.5).astype(np.int64)
    return np.atleast_2d(M)


# Search for edo maps (GPVs) that are consistent with T up to some limit
def find_edos(T, subgroup):
    assert T.ndim == 2
    r = T.shape[0]
    if r == 1:
        return

    c = kernel(T)

    # octave_div = T[0, 0]
    # print("octave mult:", octave_div)
    search_range = (4.5, 1999.5)

    m_list = []

    seen = set()
    count = 0
    count2 = 0
    for m1 in Pmaps(search_range, subgroup):
        count2 += 1
        if count2 > 20000:
            break

        # if it tempers out all commas
        if not np.any(m1 @ c):
            # if it is not contorted
            if np.gcd.reduce(m1.flatten().tolist()) == 1:
                badness = temp_badness((m1, subgroup))

                m_list.append((np.copy(m1), badness))

                # only count distinct octave divisions
                if m1[0][0] not in seen:
                    seen.add(m1[0][0])
                    count += 1
                    if count > r + 50:  # rank + 50 should be enough
                        break

    # print("list count: ", len(m_list))
    print("nr edos checked: ", count2)

    # sort by badness
    m_list.sort(key=lambda l: l[1])

    # filter so each edo only shows up once (first on the list)
    r_list = []
    seen = set()
    for m in m_list:
        if m[0][0][0] not in seen:
            r_list.append(m)
            seen.add(m[0][0][0])

    # return top 12+rank edos
    return r_list[: (r + 12)]


# Select <rank> edos that, when joined together, are equivalent to the temperament
def find_join(T, subgroup, m_list):
    assert T.ndim == 2
    r = T.shape[0]

    count = 0
    for combo in comboBySum(r, 0, len(m_list) - 1):
        m_new = np.vstack([m_list[i][0] for i in combo])
        m_hnf = hnf(m_new)

        count += 1

        if np.all(m_hnf == T):
            print("number of combos checked: " + str(count))
            return list(m_new)

        if count > 500:
            break
    print("Join search failed! Number of combos checked: " + str(count))


# Iterator for general edo maps (GPVs)
class Pmaps:
    def __init__(self, bounds, subgroup):
        self.stop = bounds[1]
        self.log_s = log_subgroup(subgroup)
        self.log_s = self.log_s / self.log_s[0]
        assert np.all(self.log_s >= 0)

        start = bounds[0]

        self.cmap = patent_map(start, subgroup)

        self.ubounds = (self.cmap + 0.5) / self.log_s

        self.first = True

    def __iter__(self):
        return self

    def __next__(self):
        if not self.first:
            incr = np.argmin(self.ubounds)
            self.cmap[0, incr] += 1
            # self.ubounds[0, incr] = (self.cmap[0, incr] + 0.5) / self.log_s[incr]
            self.ubounds[0, incr] += 1 / self.log_s[incr]

        self.first = False

        # self.lbounds = (self.cmap - 0.5) / self.log_s

        # lb = np.max(self.lbounds)
        ub = np.min(self.ubounds)

        # stop when new lower bound hits end of interval
        if ub >= self.stop:
            raise StopIteration

        # ind = (lb + ub) / 2.0
        # assert np.all(self.cmap == np.round(ind * self.log_s).astype(np.int64))

        return self.cmap


# calculates the size of some exterior product, wrt weight matrix w
# sqrt determinant of the gramian matrix
def height(M, W):
    V = M @ W
    return np.sqrt(np.linalg.det(V @ V.T))


# "logflat badness"
# Combines complexity and error into a single measurement
#
# https://en.xen.wiki/w/Tenney-Euclidean_temperament_measures#TE_logflat_badness
# The reason for the exponent here can be found in:
#
# On transfer inequalities in Diophantine approximation, II
# Y. Bugeaud, M. Laurent
# Mathematische Zeitschrift volume 265, pages249–262 (2010)
#
# See corollary 2.
# Their way of counting dimensions is different:
#  rank = d + 1
#  dim = n + 1
# Then we have
#  omega = rank / (dim - rank)
# | y ∧ X | * | X | ^ omega <= C
#
# | X | ~ complexity
# | y ∧ X | ~ error * complexity
def temp_badness(temp):
    M, S = temp
    r, d = M.shape

    j = log_subgroup(S)[np.newaxis, :]
    W = np.diagflat(1.0 / j)
    # W = np.diagflat(1.0 / np.array(S).astype(np.double))
    # W = np.eye(d)

    # normalize W so it has det(W) = 1
    W_det = np.linalg.det(W)
    assert W_det > 0
    W = W / np.power(W_det, 1 / W.shape[0])

    # the exponent
    omega = r / (d - r)

    # || M ^ j ||
    b = height(np.vstack([M, j]), W)
    # || M || ^ omega
    c = height(M, W) ** omega

    # (|| M ^ j || * || M || ^ omega ) / ||j||
    return b * c / height(j, W)
