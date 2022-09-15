# fokker block

from util import *
import numpy as np
import itertools


s = p_limit(5)
# s = [2,3,5,11]

edolist = [7,12,3]
# edolist = [9,3,2]
P = 7
# P = int(np.round(edolist[0]))


# V = np.array([[-3,2,0],[4,-1,-1],[-3,-1,2]],dtype = np.int64).T

# print()

# exit()

G = metric_weil(s)
# G= metric_tenney(s)
# G = metric_farey(7,s)
# Gd = np.linalg.inv(G)

j = log_subgroup(s)


M = np.vstack([patent_map(x,s) for x in edolist])
# M = hnf(V, transformation = True)[1]

### !!!
# M[1,3] = 9
# print(M)

# exit()

print(M)

H, U = hnf(M, transformation = True)

sol = j @ G @ H.T @ np.linalg.inv(H @ G @ H.T)
# sol, e = cte((H,s))
sol = np.atleast_2d(sol)
print(1200*sol)

print(1200*np.atleast_2d(sol) @ U)

print(M)
print(H)

print(U.T)

L = np.atleast_2d(np.array(np.round(edolist), dtype = np.int64)).T

vlist = []
for i in range(0,P+1):
	v = np.floor(0.5+i*L/P)
	# print(v)
	v = v.astype(np.int64)
	t = U@v
	tc = 1200 * sol @ t
	tc = tc.item()
	print(t.T,v.T)
	# print()
	vlist.append(tc)


for v in vlist:
	print(v)



exit()
v = M @ factors("3/2", s)
print(v)
print(H @ factors("3/2", s))
print(U @ v)