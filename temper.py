# collection of functions for dealing with regular temperaments

import olll
from fractions import Fraction
import numpy as np
import diophantine
from math import gcd
from primes import primes
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

# expansion for fractional subgroup
def prime_expansion(subgroup):
	assert (type(subgroup) is list), "Subgroup must be a list."
	assert all(isinstance(x, Fraction) for x in subgroup)

	s = set()

	for p in subgroup:
		s = s.union(prime_factors(p))

	return sorted(list(s))


def prime_factors(frac):
	assert type(frac) is Fraction

	s = set()

	n,d = frac.as_integer_ratio()

	for p in primes:
		# print(p,n,d)
		while n % p == 0:
			# print(p)
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


def factors(fr, subgroup):
	assert (type(subgroup) is list), "Subgroup must be a list."
	assert all(isinstance(x, int) for x in subgroup)

	factors = np.zeros((len(subgroup), 1), dtype=np.int64)

	if type(fr) is tuple:
		p = fr[0]
		q = fr[1]
	else:
		frac = Fraction(fr)
		p = frac.numerator
		q = frac.denominator

	assert (p > 0)
	assert (q > 0)

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


def make_positive(v, subgroup):
	if log_interval(v, subgroup) < 0:
		return -v
	else:
		return v

def log_subgroup(subgroup):
	return np.log2(np.array(subgroup).astype(np.double))

def log_interval(v, subgroup):
	vn = v.flatten().tolist()
	logs = log_subgroup(subgroup)
	r = 0
	for i, n in enumerate(vn):
		r += n * logs[i]

	return r


def ratio(v, subgroup):
	# assert all(isinstance(x, Fraction) for x in subgroup)

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
	logs = log_subgroup(subgroup)
	# floor(x+0.5) rounds more predictably (downwards on .5)
	M = np.floor(t * logs + 0.5).astype(np.int64)
	return np.atleast_2d(M)


def p_limit(limit):
	assert limit > 1
	assert limit < 7920

	return [i for i in primes if i <= limit]


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
			# print(M)
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

	M = np.atleast_2d(M)

	j = log_subgroup(s)

	W = np.diag(1/j)
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

	W = np.diag(1/j)
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


# patent map iterator
class Pmaps:
	def __init__(self, bounds, subgroup):
		self.stop = bounds[1]
		self.logS = log_subgroup(subgroup)
		assert np.all(self.logS >= 1)

		start = bounds[0]

		self.cmap = patent_map(start, subgroup)

		self.first = True

	def __iter__(self):
		return self

	def __next__(self):
		if not self.first:
			incr = np.argmin(self.ubounds)
			self.cmap[0, incr] += 1

		self.first = False

		self.lbounds = (self.cmap - 0.5) / self.logS
		self.ubounds = (self.cmap + 0.5) / self.logS

		lb = np.max(self.lbounds)
		ub = np.min(self.ubounds)

		# stop when new lower bound hits end of interval
		if lb >= self.stop:
			raise StopIteration

		ind = (lb + ub) / 2.0
		assert np.all(self.cmap == np.round(ind * self.logS).astype(np.int64))

		return self.cmap, (lb, ub)


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

def get_subgroup_basis(subgroup):
	expanded = prime_expansion(subgroup)
	s_basis = []
	for k in subgroup:
		s_basis.append(factors(k,expanded))
	s_basis = np.hstack(s_basis).T
	
	# new_subgroup = subgroup

	# if redundant, normalize
	if np.linalg.det(s_basis @ s_basis.T) < 0.5:
		s_basis = hnf(s_basis, remove_zeros = True)
		# new_subgroup = []
		# for b in s_basis:
			# new_subgroup.append(ratio(b,expanded))

	return s_basis, expanded

def get_subgroup(s_basis, expanded):
	s = []
	for b in s_basis:
		s.append(ratio(b,expanded))
	return s

if __name__ == '__main__':
	subgroup = [Fraction("2"),Fraction("5/4"),Fraction("9")]

	s_basis, expanded = get_subgroup_basis(subgroup)

	subgroup = get_subgroup(s_basis, expanded)

	rational = False
	for r in subgroup:
		p,q = r.as_integer_ratio()
		if q > 1:
			rational = True

	if not rational:
		subgroup = list(map(int, subgroup))


	# normalization test.. unsuited because it doesnt keep equave
	# s_basis = hnf(s_basis, remove_zeros = True)
	# s_basis = LLL(s_basis)
	# # make positive
	# for i in range(s_basis.shape[0]):
	# 	if log_interval(s_basis[i], expanded) < 0:
	# 		s_basis[i, :] = -s_basis[i, :]
	# new_subgroup = []
	# for b in s_basis:
	# 	new_subgroup.append(ratio(b,expanded))

	print(expanded)
	print(s_basis)
	print(subgroup)
