# https://github.com/higumachan/olll
# corrected, and modified to use an arbitrary positive definite bilinear form

import numpy as np

from .util_types import IntMat, FloatMat, IntVec

MAX_LLL_ITERS = 300


def innerprod(a: IntVec, b: IntVec, W: FloatMat) -> float:
    return np.dot(a, W @ b)


def gramschmidt(v: IntMat, W: FloatMat) -> FloatMat:
    v = v.astype(np.double)
    u = v.copy()

    for i in range(1, v.shape[0]):
        ui = u[i]
        for uj in u[:i]:
            ui = ui - (innerprod(uj, v[i], W) / innerprod(uj, uj, W)) * uj

        u[i] = ui
    return u


# LLL reduction on the rows of `basis`
def reduction(basis: IntMat, delta: float, W: FloatMat) -> IntMat:
    n = basis.shape[0]
    ortho = gramschmidt(basis, W)

    def mu(i: int, j: int) -> float:
        a, b = ortho[j], basis[i]
        return innerprod(a, b, W) / innerprod(a, a, W)

    iter_count = 0

    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            mu_kj = mu(k, j)
            if abs(mu_kj) > 0.5:
                basis[k] = basis[k] - basis[j] * round(mu_kj)
                ortho = gramschmidt(basis, W)

        l_condition = (delta - mu(k, k - 1) ** 2) * innerprod(ortho[k - 1], ortho[k - 1], W)
        if innerprod(ortho[k], ortho[k], W) >= l_condition:
            k += 1
        else:
            # basis[k], basis[k - 1] = basis[k - 1], basis[k].copy()  ##fix
            basis[[k - 1, k]] = basis[[k, k - 1]]
            ortho = gramschmidt(basis, W)
            k = max(k - 1, 1)

        if iter_count > MAX_LLL_ITERS:
            # we are likely stuck in an infinite loop due to precision issues
            raise OverflowError("LLL reduction took too long")
        iter_count += 1
    return basis


# Babai's nearest plane algorithm for solving approximate CVP
# `basis` should be LLL reduced first
# operates on rows
def nearest_plane(v: IntVec, basis: IntMat, W: FloatMat) -> IntVec:
    b = v.copy()
    n = basis.shape[0]

    ortho = gramschmidt(basis, W)

    for j in range(n - 1, -1, -1):
        a = ortho[j]
        mu = innerprod(a, b, W) / innerprod(a, a, W)
        b = b - round(mu) * basis[j]
    return v - b
