# collection of functions for dealing with regular temperaments

from . import olll
from fractions import Fraction
import numpy as np
# from math import gcd
# from primes import primes
from itertools import combinations
from . import diophantine

from .subgroup import *
from .optimize import *
from .combo import comboBySum



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


def LLL(M, W, delta = 0.99):
	res = olll.reduction(np.copy(M).T, delta=delta, W=W).T
	# res2 = np.array(olll2.reduction(np.copy(M).T, delta=0.99), dtype=np.int64).T

	# sort 'em
	# actually, this might be redundant.
	# c_list = list(res.T)
	# c_list.sort(key = lambda c: np.dot(c, W @ c))
	# return np.array(c_list).T
	return res


def antitranspose(M):
	return np.flipud(np.fliplr((M.T)))


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


def factor_order(M):
	r, d = M.shape
	return integer_det(hnf(M.T)[:r].T)


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
	assert np.all(A @ sol == B), "Could not solve system"

	return sol


def preimage(M):
	gens = []

	r, d = M.shape

	rank = M.shape[0]

	gens = solve_diophantine(M, np.eye(rank, dtype=np.int64))

	return gens

# simplify Intervals wrt Comma basis (should be in reduced LLL) wrt Weight matrix
def simplify(I,C,W):
	It = I.T
	Ct = C.T

	for i in range(len(It)):
		v = It[i]
		p_best = np.dot(v, W @ v)

		cont = True
		while cont:
			# print(v, np.dot(v, W @ v))
			cont = False
			for c in Ct:
				new = v - c
				p_new = np.dot(new, W @ new)
				if p_new < p_best:
					v = new
					p_best = p_new
					cont = True
				else:
					new = v + c
					p_new = np.dot(new, W @ new)
					if p_new < p_best:
						v = new
						p_best = p_new
						cont = True
		It[i] = v

	return It.T

def patent_map(t, subgroup):
	logs = log_subgroup(subgroup)

	t = t / logs[0]  # fix equave

	# floor(x+0.5) rounds more predictably (downwards on .5)
	M = np.floor(t * logs + 0.5).astype(np.int64)
	return np.atleast_2d(M)

def find_edos_patent(T, subgroup):
	assert (T.ndim == 2)
	r, d = T.shape
	# T = hnf(T)
	c = kernel(T)

	octave_div = T[0,0]
	# print("octave mult:", octave_div)
	# search_range = (4.5, 665.5)

	m_list = []

	if r == 1:
		return

	seen = set()
	count = 0
	count2 = 0
	for k in range(666):
		m1 = patent_map(k,subgroup)
		# print(m1[0,0])
		if m1[0,0] % octave_div == 0: # skip non multiples of the octave division
			count2 += 1
			if count2 > 8000:
				break
			# if it tempers out all commas
			if np.all(m1 @ c == 0):
				# if it is not contorted
				if np.gcd.reduce(m1.flatten().tolist()) == 1:
					badness = temp_measures((m1, subgroup))[0]
					m_list.append((np.copy(m1), badness))

					# only count distinct octave divisions
					if m1[0][0] not in seen:
						seen.add(m1[0][0])
						count += 1
						if count > r + 10:  # rank + 10 should be enough
							break


	print("list count: ", len(m_list))
	print("nr checked: ", count2)

	# sort by badness
	m_list.sort(key=lambda l: l[1])

	# filter so each edo only shows up once (first on the list)
	r_list = []
	seen = set()
	for m in m_list:
		if m[0][0][0] not in seen:
			r_list.append(m)
			seen.add(m[0][0][0])

	return r_list

def find_edos(T, subgroup):
	assert (T.ndim == 2)
	r, d = T.shape
	# T = hnf(T)
	c = kernel(T)

	octave_div = T[0,0]
	# print("octave mult:", octave_div)
	search_range = (4.5, 1999.5)

	m_list = []

	if r == 1:
		return

	seen = set()
	count = 0
	count2 = 0
	for m1, b1 in Pmaps(search_range, subgroup):
		# print(m1[0,0])
		if m1[0,0] % octave_div == 0: # skip non multiples of the octave division
			count2 += 1
			if count2 > 8000:
				break
			# if it tempers out all commas
			if np.all(m1 @ c == 0):
				# if it is not contorted
				if np.gcd.reduce(m1.flatten().tolist()) == 1:
					badness = temp_measures((m1, subgroup))[0]

					patent = (np.floor(b1[0]) == np.floor(b1[1]))
					m_list.append((np.copy(m1), badness, patent))

					# only count distinct octave divisions
					if m1[0][0] not in seen:
						seen.add(m1[0][0])
						count += 1
						if count > r + 25:  # rank + 10 should be enough
							break


	print("list count: ", len(m_list))
	print("nr checked: ", count2)

	# sort by badness
	m_list.sort(key=lambda l: l[1])

	# filter so each edo only shows up once (first on the list)
	r_list = []
	seen = set()
	for m in m_list:
		if m[0][0][0] not in seen:
			r_list.append(m)
			seen.add(m[0][0][0])

	# return r_list
	return r_list[:(r+12)]



def find_join(T, subgroup, m_list):
	assert (T.ndim == 2)
	r, d = T.shape
	# T = hnf(T)

	count = 0
	for combo in comboBySum(r, 0, len(m_list) - 1):
		# print(combo, flush=True)
		m_new = np.vstack([m_list[i][0] for i in combo])
		m_hnf = hnf(m_new)

		# print(m_hnf)
		count += 1

		if (np.all(m_hnf == T)):
			print("number of combos checked: " + str(count))
			return [m[0] for m in m_new]

		if count > 500:
			break
	print("FAILED. number of combos checked: " + str(count))

# patent map iterator
class Pmaps:
	def __init__(self, bounds, subgroup):
		self.stop = bounds[1]
		self.logS = log_subgroup(subgroup)
		# assert np.all(self.logS >= 1)

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


def temp_relerror(temp, G):
	M, S = temp
	r, d = M.shape

	j = np.atleast_2d(log_subgroup(S))

	# js = np.sqrt(np.asscalar(j @ G @ j.T))

	V = np.block([[M],[j]])


	err = np.sqrt(np.linalg.det(V @ G @ V.T))

	# sol, e = lstsq(temp, weight="tenney")
	# err2 = np.sqrt(np.average(e**2))
	# print(np.sqrt((d-r)/((r+1)*d))*err2)

	# Breed
	# err2 = np.sqrt(np.average((e @ W)**2))


	# Smith
	# err2 = np.sqrt(np.sum((e @ W)**2) * (r + 1) / (d-r) )

	return err


def temp_complexity(temp, G):
	M, S = temp
	r, d = M.shape


	compl = np.sqrt(np.linalg.det((M @ G @ M.T)))

	# simple
	# compl = np.sqrt(np.linalg.det((M @ G @ M.T)) / d)

	# Breed
	# compl  = np.sqrt (np.linalg.det ((M @ W) @ (M @ W).T / d))

	# Smith
	# compl = np.sqrt (np.linalg.det ((M @ W) @ (M @ W).T ) / math.comb(d,r))

	return compl


# logflat badness
epsilon = 0.00
def temp_measures(temp, G):
	M, S = temp
	r, d = M.shape

	complexity = temp_complexity(temp, G)

	relerror = temp_relerror(temp, G)
	error = relerror / complexity

	badness = error * (complexity**(epsilon + (d / (d - r))))

	# c2 = (complexity**(epsilon + (d / (d - r))))

	# print(relerror, complexity**(- r / (d - r)))

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
