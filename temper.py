# collection of functions for dealing with regular temperaments

import olll
from fractions import Fraction
import numpy as np
import diophantine
from math import gcd

# PHI = (1 + 5 ** 0.5) / 2


def hnf(M, remove_zeros=False):
	assert (M.ndim == 2)

	res = np.array(diophantine.lllhermite(M.astype(np.int64))[0]).astype(np.int64)

	if remove_zeros:
		idx = np.argwhere(np.all(res[:] == 0, axis=1))
		res = np.delete(res, idx, axis=0)

	return res


def kernel(M):
	assert (M.ndim == 2)
	r, d = M.shape
	n = d - r

	M = (np.vstack([M, np.eye(d, dtype=np.int64)]))
	K = hnf(M.T).T[r::, r::]

	return K


def cokernel(M):
	assert (M.ndim == 2)
	return kernel(M.T).T


def LLL(M):
	return np.array(olll.reduction(M.T.tolist(), delta=0.99), dtype=np.int64).T


def antitranspose(M):
	return np.flipud(np.fliplr((M.T)))


def defactored_hnf(M):
	r, d = M.shape

	S = np.linalg.inv(hnf(M.T)[:r].T)

	assert np.allclose(S @ M, np.round(S @ M))

	D = np.round(S @ M).astype(np.int64)

	return hnf(D)


def factor_order(M):
	r, d = M.shape
	return np.round(np.linalg.det(hnf(M.T)[:r].T)).astype(np.int64)


def canonical(M):
	assert (M.ndim == 2)
	r, d = M.shape
	if r > d:
		# comma basis
		return antitranspose(defactored_hnf(antitranspose(M)))
	else:
		# mapping
		return defactored_hnf(M)


def prime_factors(n, subgroup):
	factors = np.zeros((len(subgroup), 1), dtype=np.int64)

	for i, p in enumerate(subgroup):
		while n % p == 0:
			n //= p
			factors[i, 0] += 1

	assert n == 1, "Decomposition not in subgroup."

	return factors


def factors(fr, subgroup):
	assert (type(subgroup) is list), "Subgroup must be a list."
	if type(fr) is tuple:
		p = fr[0]
		q = fr[1]
	else:
		frac = Fraction(fr)
		p = frac.numerator
		q = frac.denominator

	assert (p >= 0)
	assert (q > 0)

	return prime_factors(p, subgroup) - prime_factors(q, subgroup)


def make_positive(v, subgroup):
	if log_interval(v, subgroup) < 0:
		return -v
	else:
		return v


def log_interval(v, subgroup):
	vn = v.flatten().tolist()
	r = 0
	for i, n in enumerate(vn):
		r += n * np.log2(subgroup[i])

	return r


def ratio(v, subgroup):
	v = v.flatten().tolist()
	p = 1
	q = 1
	for i, n in enumerate(v):
		if n > 0:
			p = p * subgroup[i]**n
		elif n < 0:
			q = q * subgroup[i]**(-n)

	return Fraction(p, q)


def patent_map(t, subgroup):
	factors = np.zeros((1, len(subgroup)), dtype=np.int64)

	M = np.zeros((1, len(subgroup)), dtype=np.int64)

	for i, n in enumerate(subgroup):
		M[0, i] = round(t * np.log2(n))

	return M


def p_limit(limit):
	limitn = limit + 1
	primes = dict()
	for i in range(2, limitn):
		primes[i] = True

	for i in primes:
		factors = range(i, limitn, i)
		for f in factors[1:]:
			primes[f] = False
	return [i for i in primes if primes[i] == True]


def findMaps(T, subgroup):
	assert (T.ndim == 2)
	r, d = T.shape

	c = 0
	Tp = np.zeros((1, d), dtype=np.int64)

	f_order = factor_order(T)

	for i in np.arange(5, 500, 1):
		M = (patent_map(i, subgroup))

		# check if edo map supports given M
		red = hnf(np.vstack([T, M]))[r, :]

		# check if edo map has gcd = 1
		tors = np.gcd.reduce(M.flatten().tolist())

		if not red.any() and tors == 1:
			# check if new edo map is not redundant
			newmap = hnf(np.vstack([Tp, M]))
			red2 = newmap[c, :]

			# check if contorsion matches input
			order = factor_order(newmap)

			if (red2.any() and order == f_order) or c == 0:
				if c == 0:
					Tp = M
				else:
					Tp = np.vstack([Tp, M])

				# print(Tp)
				c = c + 1
				if c >= r:
					break
	return Tp


def preimage(M):
	gens = []

	r, d = M.shape

	for i in range(r):
		tv = np.zeros(r, dtype=np.int64)
		tv[i] = 1
		sol = diophantine.solve(M, tv)
		g = sol[0]
		gens.append(g)

	return gens


def lstsq(temp, weight="tenney"):
	M = temp[0]
	s = temp[1]

	print(M)

	j = np.log2(np.array(s))

	W = np.diag(1. / j)
	if weight == 'unweighted':
		W = np.eye(len(s))

	j = np.atleast_2d(j)

	sol = np.linalg.lstsq((M @ W).T, (j @ W).T, rcond=None)[0]

	tun = (sol.T @ M)
	err = tun - j

	return sol.flatten(), err.flatten()


def cte(temp, weight="tenney"):
	M = temp[0]
	s = temp[1]

	j = np.log2(np.array(s))

	W = np.diag(1. / j)
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


# def normsq(v,W = None):
#     if W is None:
#         return np.sum(v*v)
#     else:
#         return np.sum(v*W*v)

# def norm(v,W = None):
#     return np.sqrt(normsq(v,W))

# def inner(a,b,W = None):
#     if W is None:
#         return np.sum(a*b)
#     else:
#         return np.sum(a*W*b)
