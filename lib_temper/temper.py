# collection of functions for dealing with regular temperaments

from typing import Optional, TypeVar

import numpy as np
import rslattice

from .util_types import IntMat, FloatMat
from .subgroup import *
from .optimize import *
from .combo import comboBySum


# Find the hermite normal form of M
def hnf(M: IntMat, remove_zeros: bool = False) -> IntMat:
    assert M.ndim == 2

    res = rslattice.hnf(M)

    if remove_zeros:
        idx = np.argwhere(np.all(res[:] == 0, axis=1))
        res = np.delete(res, idx, axis=0)

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


# LLL reduction
# Tries to find the smallest basis vectors in some lattice
# Operates on columns
# M: Map
# W: Weight matrix (metric)
def LLL(M: IntMat, W: Optional[FloatMat] = None, delta: float = 0.75) -> IntMat:
    if W is None:
        W = np.eye(M.shape[0])

    res = rslattice.lll(M.T, delta, W).T

    # sort them by complexity
    # actually, this might be redundant.
    c_list = list(res.T)
    c_list.sort(key=lambda c: np.dot(c, W @ c))
    return np.array(c_list).T  # type: ignore[return-value]

    # return res


# Flips M along diagonal
def antitranspose(M: IntMat) -> IntMat:
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
def defactored_hnf(M: IntMat) -> IntMat:
    r, _ = M.shape

    K = rslattice.hnf(M.T)[:r].T
    if rslattice.integer_det(K) == 1:
        return rslattice.hnf(M)

    S = np.linalg.inv(K)

    D = np.round(S @ M)
    assert np.allclose(S @ M, D)

    return hnf(D.astype(np.int64))


# Order of factorization.
# For a saturated basis this is 1.
def factor_order(M: IntMat) -> int:
    r = M.shape[0]
    return rslattice.integer_det(rslattice.hnf(M.T)[:r].T)


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
def solve_diophantine(A, B):
    assert A.shape[0] == B.shape[0]
    _, d = A.shape

    B = np.atleast_2d(B)
    aug = np.block([[A, B]])

    H = rslattice.hnf(aug)
    sol = H[0:d, d:]

    K = H[0:d, 0:d]
    S = np.linalg.inv(K)
    sol = np.round(S @ sol).astype(np.int64)

    # Check that the solution actually works.
    # Probably the easiest way to guarantee this routine works correctly.
    assert np.all(A @ sol == B), "Could not solve system"

    return sol


# Solves XM = I in the integers (basically as above)
def preimage(M):
    r, d = M.shape

    B = np.block([[M.T, np.eye(d, dtype=np.int64)]])
    H = rslattice.hnf(B)
    sol = H[0:r, r:].T

    assert np.all(M @ sol == np.eye(r, dtype=np.int64)), "Could not solve system"

    return sol


# Simplify a list of intervals with respect to some comma basis.
# The comma basis should be in reduced LLL form for this to work properly.
# To find the solution, we solve approximate CVP in the lattice spanned by the commas.
def simplify(intervals, commas, W=None):
    if W is None:
        W = np.eye(commas.shape[0])

    simpl = intervals.copy()
    for i in range(simpl.shape[1]):
        simpl[:, i] -= rslattice.nearest_plane(simpl[:, i], commas.T, W)

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
    g = np.linalg.det(M @ W @ M.T)
    if g <= 1e-8:
        return 1e-4

    return np.sqrt(g)


# "logflat badness" aka Dirichlet coefficient
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
def temp_badness(temp, W=None):
    M, S = temp
    r, d = M.shape

    j = log_subgroup(S)[np.newaxis, :]
    if W is None:
        W = np.diagflat(1.0 / j) ** 2

    # W = np.eye(d)
    # normalize W so it has det(W) = 1
    W_det = np.linalg.det(W)
    assert W_det > 0
    W = W / np.power(W_det, 1 / d)

    # W = W / np.power(W_det, 1 / r)
    # W_det = np.linalg.det(W)

    # the exponent
    omega = r / (d - r)

    # || M ^ j ||
    b = height(np.vstack([M, j]), W)
    # b = height(np.vstack([M, j]), W) / ((W_det) ** (1/2))
    # || M || ^ omega
    c = height(M, W)
    # c = height(M, W) / ((W_det) ** (1/2))
    c = c**omega


    # (|| M ^ j || * || M || ^ omega ) / ||j||
    return b * c / height(j, W)
    # return b * c
    # return c
