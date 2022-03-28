## generator optimization

import numpy as np

from .interval import *


def lstsq(temp, weight="tenney"):

	M = temp[0]
	s = temp[1]

	M = np.atleast_2d(M)

	j = log_subgroup(s)

	W = np.diag(1 / j)
	if weight == 'unweighted':
		W = np.eye(len(s))

	j = np.atleast_2d(j)

	sol = np.linalg.lstsq((M @ W).T, (j @ W).T, rcond=None)[0]

	tun = (sol.T @ M)
	err = tun - j

	return sol, err.flatten()


def cte(temp, weight="tenney"):
	M = temp[0]
	s = temp[1]
	M = np.atleast_2d(M)

	j = log_subgroup(s)

	W = np.diag(1 / j)
	if weight == 'unweighted':
		W = np.eye(len(s))

	j = np.atleast_2d(j)

	A = (M @ W).T
	b = (j @ W).T

	r, d = M.shape

	V = np.zeros((d, 1))
	V[0, 0] = 1

	C = (M @ V).T
	d = (j @ V).T

	sol = np.linalg.solve(np.block([[A.T @ A, C.T], [C, 0]]), np.vstack([A.T @ b, d]))

	sol = sol[:r]

	tun = (sol.T @ M)
	err = tun - j

	return sol.flatten(), err.flatten()