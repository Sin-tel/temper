# https://github.com/higumachan/olll
# corrected, and modified to use an arbitrary positive definite bilinear form

from typing import Sequence, List

import numpy as np


def innerprod(a, b, W):
	return np.dot(a, W @ b)


def gramschmidt(v: np.ndarray, W) -> np.ndarray:
	v = v.astype(np.double)
	u = v.copy()

	for i in range(1, len(v)):
		ui = u[i]
		for uj in u[:i]:
			ui = ui - (innerprod(uj, v[i], W) / innerprod(uj, uj, W)) * uj

		u[i] = ui
	return u


def reduction(basis: np.ndarray, delta: float, W) -> Sequence[Sequence[int]]:
	n = len(basis)
	ortho = gramschmidt(basis, W)

	def mu(i: int, j: int) -> float:
		a, b = ortho[j], basis[i]
		return innerprod(a, b, W) / innerprod(a, a, W)

	k = 1
	while k < n:
		for j in range(k - 1, -1, -1):
			mu_kj = mu(k, j)
			if abs(mu_kj) > 0.5:
				basis[k] = basis[k] - basis[j] * round(mu_kj)
				ortho = gramschmidt(basis, W)

		if innerprod(ortho[k], ortho[k],
					 W) >= (delta - mu(k, k - 1)**2) * innerprod(ortho[k - 1], ortho[k - 1], W):
			k += 1
		else:
			basis[k], basis[k - 1] = basis[k - 1], basis[k].copy()  ##fix
			ortho = gramschmidt(basis, W)
			k = max(k - 1, 1)

	return basis
