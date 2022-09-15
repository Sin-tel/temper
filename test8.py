from util import *
import numpy as np
import math

s = p_limit(5)

j = log_subgroup(s)
W = np.diag(1.0 / j)
# W[0,0] = 0

print(W)
def metric_kees(s):
	s = s[1:]
	n = len(s)
	j = log_subgroup(s)
	Bw = np.block([[np.eye(n)],[np.ones((1,n))]])
	# Bw = np.eye(n)

	Bw2 = Bw @ np.diag(j.flatten())
	# print("===")
	# print(Bw2)
	Gd = Bw2.T @ Bw2

	G = np.linalg.inv(Gd)
	# G = G / G[0,0]

	z = np.zeros((len(s),1))
	G = np.block([[0,z.T],[z,G]])
	Gd = np.block([[0,z.T],[z,Gd]])

	# print(Gd)

	return G, Gd

def cte2(temp, weight="tenney", V = None):
	M = temp[0]
	s = temp[1]
	M = np.atleast_2d(M)

	j = log_subgroup(s)
	j = np.atleast_2d(j)

	r, d = M.shape

	if V is None:
		V = np.zeros((d, 1))
		V[0, 0] = 1

	C = (M @ V).T
	d = (j @ V).T

	# remove octave
	# M = M[:,1:]
	# j = j[:,1:]


	A = (M).T
	b = (j).T




	Z = np.zeros((C.shape[0], C.shape[0]))
	# G, Gd = metric_kees(s)
	# G = G[1:,1:]
	# G = metric_weil(s[1:])
	G = metric_weil(s)
	
	# W = np.diag(1/j.flatten())
	# G = W @ W
	# print(metric_weil(s[1:]))
	# print(G/metric_weil(s[1:]))
	

	sol = np.linalg.solve(np.block([[A.T @ G @ A, C.T], [C, Z]]), np.vstack([A.T @ G @ b, d]))

	sol = sol[:r]

	tun = (sol.T @ M)
	err = tun - j

	return sol.flatten(), err.flatten()

# M = 
e12 = patent_map(12,s)
e19 = patent_map(19,s)
M = np.vstack([e12,e19])
M = hnf(M)

print(M)

j = np.atleast_2d(j)

G = W @ W.T
sol = j @ G @ M.T @ np.linalg.inv(M @ G @ M.T)
sol = sol.flatten()
sol /= sol[0]*M[0,0]

print(1200*sol)

sol2, e = cte((M,s), weight="tenney")
sol2 = sol2.flatten()
# sol2 /= sol2[0]
print(1200*sol2)

sol3, e = cte2((M,s), weight="tenney")
sol3 = sol3.flatten()
# sol3 /= sol3[0]
print(1200*sol3)

# G, Gd = metric_kees(s)
G = metric_wilson(s)
sol4 = j @ G @ M.T @ np.linalg.inv(M @ G @ M.T)
sol4 = sol4.flatten()
sol4 /= sol4[0]*M[0,0]

print(1200*sol4)
# exit()
# Gd = Gd/2
print("===============================")


G, Gd = metric_kees(s)
# G = metric_tenney(s)
# Gd = np.linalg.inv(G)
tmap = sol

print(1200*tmap)

errmap = ((tmap@M - j))
# e = (tmap@M - j)
print(1200*errmap)
# exit()
# errm = np.sqrt(errmap @ G @ errmap.T)
errm = (errmap/j).flatten()[1:]
# exit()
errm = (np.max(errm) - np.min(errm))
# print(e)
# exit()
# errm = np.max(errmap.flatten()[1:])
# errm = errmap/j
print(errm)
# exit()
print("norm of error map:")
print(1200*errm)
print("================")
# exit()

# v1 = factors("2/1",s)
# e1 = 1200*(tmap @ M @ v1 - log_interval(v1,s))

print(1200 * sol @ M)
print(1200 * sol3 @ M)
# print(e1)
exit()

l = farey(1000)
lv = set()
for i in l:
	f = factors_unchecked(i, s)
	if f is not None:
		f[0,0] = 0
		if not np.all(f == 0):
			lv.add(tuple(f.flatten().tolist()))

total = 0
maxe = 0
lv  = list(lv)
for k in lv:
	i = np.array(k)
	i = np.atleast_2d(i).T

	frac = ratio(i,s)
	p = frac.numerator
	q = frac.denominator
	# w = np.log2(max(p,q))

	# w = np.sqrt((i.T @ Gd @ i))

	print(ratio(i,s))
	w = (i.T*j).flatten()[1:]
	print(w)
	# print(i*j.T)
	w = 0.5 * (np.sum(np.abs(w)) + np.abs(np.sum(w)))
	print(w)
	print(np.log2(max(p,q)))

	# print(np.log2(max(p,q)))
	# print(w)
	# w = np.sum(np.abs(w))
	if w > 0:

		e = 1200*(tmap @ M @ i - log_interval(i,s))

		e_adj = np.abs(e)/w
		maxe = max(maxe,e_adj) 
		print(e_adj)
		total += e_adj


print(len(lv))
print("maxerror:")
print(maxe)
print("error:")
print(total / len(lv))
print("norm of error map:")
print(1200*errm)
# 8.07161105