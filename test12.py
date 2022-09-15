from util import *
import numpy as np
import math

s = p_limit(7)

j = log_subgroup(s)


e12 = patent_map(12,s)
e19 = patent_map(19,s)
M = np.vstack([e12,e19])
M = hnf(M)

T = M
gens = preimage(T)
# eq = log_interval(gens[0], s)
o = T[0, 0]
genoct = np.zeros_like(gens[:, 0])
genoct[0] = 1
# reduce by octave
for i in range(1, T.shape[0]):
	# make positive first
	if log_interval(gens[:, i], s) < 0:
		T[i, :] = -T[i, :]
		gens[:, i] = -gens[:, i]

	red = int(np.floor(log_interval(gens[:, i], s)))
	gens[:, i] -= red * genoct
	T[0, :] += o * red * T[i, :]

M = T
# M[0,2] = 10
print(M)


G = metric_weil(s)
# G = metric_farey(6,s)

sol = j @ G @ M.T @ np.linalg.inv(M @ G @ M.T)
# sol = sol.flatten()


print(1200*sol)


print(1200*np.atleast_2d(sol) @ M)
print(1200*(np.atleast_2d(sol) @ M - j))