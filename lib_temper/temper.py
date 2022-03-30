# collection of functions for dealing with regular temperaments

import olll
from fractions import Fraction
# import numpy as np
# from math import gcd
# from primes import primes

from . import diophantine

from .subgroup import *
from .optimize import *


def hnf(M, remove_zeros=False, transformation=False):
	assert (M.ndim == 2)

	solution = diophantine.lllhermite(M.astype(np.int64))

	res = np.array(solution[0]).astype(np.int64)

	unimod = np.array(solution[1]).astype(np.int64)

	if remove_zeros:
		idx = np.argwhere(np.all(res[:] == 0, axis=1))
		res = np.delete(res, idx, axis=0)

	if transformation:
		return res, unimod
	else:
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


# Solve AX = B in the integers
# for the method used, see https://github.com/tclose/Diophantine/blob/master/algorithm.pdf
def solve_diophantine(A, B):
	B = np.atleast_2d(B)
	assert A.shape[0] == B.shape[0]
	aug = np.block([[A.T, np.zeros((A.shape[1], B.shape[1]), dtype=np.int64)],
	                [B.T, np.eye(B.shape[1], dtype=np.int64)]])

	r, d = A.shape

	nullity = d - r

	## somehow we can solve the system even if the nullity is -1
	# aka the kernel is trivial
	# should double check when this actually works
	if nullity <= 0:
		nullity = 0

	H, U = hnf(aug, transformation=True)

	p2 = U.shape[0] - nullity
	p1 = p2 - B.shape[1]

	sol = -U[p1:p2, :A.shape[1]].T

	# Check that the solution actually works.
	# Probably the easiest way to guarantee this routine works correctly.
	assert np.all(A @ sol == B)

	return sol.T


def preimage(M):
	gens = []

	r, d = M.shape

	rank = M.shape[0]

	gens = solve_diophantine(M, np.eye(rank, dtype=np.int64))

	return gens


def patent_map(t, subgroup):
	logs = log_subgroup(subgroup)

	t = t / logs[0]  # fix equave

	# floor(x+0.5) rounds more predictably (downwards on .5)
	M = np.floor(t * logs + 0.5).astype(np.int64)
	return np.atleast_2d(M)


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


def temp_error(temp):
	M, S = temp
	r, d = M.shape

	j = log_subgroup(S)
	W = np.diag(1. / j)

	sol, e = lstsq(temp, weight="tenney")

	# Breed
	err = np.sqrt(np.average((e @ W)**2))

	# Smith
	# err = np.sqrt(np.sum((e @ W)**2) * (r + 1) / (d-r) )

	return err


def temp_complexity(temp):
	M, S = temp
	r, d = M.shape

	j = log_subgroup(S)
	W = np.diag(1. / j)

	# simple
	compl = np.sqrt(np.linalg.det((M @ W) @ (M @ W).T) / d)

	# Breed
	# compl  = np.sqrt (np.linalg.det ((M @ W) @ (M @ W).T / d))

	# Smith
	# compl = np.sqrt (np.linalg.det ((M @ W) @ (M @ W).T ) / math.comb(d,r))

	return compl


# logflat badness
def temp_measures(temp):
	M, S = temp
	r, d = M.shape

	complexity = temp_complexity(temp)
	error = temp_error(temp)

	badness = error * (complexity**(d / (d - r)))

	return badness, complexity, error


if __name__ == '__main__':
	subgroup = [Fraction("2"), Fraction("5/4"), Fraction("9")]

	s_basis, expanded = get_subgroup_basis(subgroup)

	subgroup = get_subgroup(s_basis, expanded)

	rational = False
	for r in subgroup:
		p, q = r.as_integer_ratio()
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
