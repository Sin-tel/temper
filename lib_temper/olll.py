# https://github.com/higumachan/olll
# corrected, and modified to use an arbitrary positive definite bilinear form

import numpy as np

from .util_types import IntMat, FloatMat, IntVec


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


def reduction(basis: IntMat, delta: float, W: FloatMat) -> IntMat:
    n = basis.shape[0]
    ortho = gramschmidt(basis, W)

    type_info = np.iinfo(np.int64)
    i64_min = type_info.min
    i64_max = type_info.max

    def mu(i: int, j: int) -> float:
        a, b = ortho[j], basis[i]
        return innerprod(a, b, W) / innerprod(a, a, W)

    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            mu_kj = mu(k, j)
            if abs(mu_kj) > 0.5:
                if mu_kj < i64_min or mu_kj > i64_max:
                    raise OverflowError("your numbers are too big!")

                basis[k] = basis[k] - basis[j] * np.round(mu_kj).astype(np.int64)
                ortho = gramschmidt(basis, W)

        l_condition = (delta - mu(k, k - 1) ** 2) * innerprod(ortho[k - 1], ortho[k - 1], W)
        if innerprod(ortho[k], ortho[k], W) >= l_condition:
            k += 1
        else:
            # basis[k], basis[k - 1] = basis[k - 1], basis[k].copy()  ##fix
            basis[[k - 1, k]] = basis[[k, k - 1]]
            ortho = gramschmidt(basis, W)
            k = max(k - 1, 1)

    return basis
