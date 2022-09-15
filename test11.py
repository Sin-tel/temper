from util import *
import numpy as np
import math
import itertools
# from time import perf_counter

s = p_limit(7)
# s = [2,3/2,5/4,7/4,11/8]
# s = [2,3,7]

n = len(s)
j = log_subgroup(s)

search_range = (5,300)

G = metric_weil(s)
# G = metric_tenney(s)
# G = np.eye(len(s))
# print(np.linalg.inv(G))

# G = G / (np.dot(j,G@j)) 
# print(1/n)
# print(np.linalg.det(G)**(1/n))
# print(np.dot(j,G@j))
# print(np.linalg.det(G))
# exit()
# G = G / (np.linalg.det(G)**(1/n)) # -> normalizes for rank n
G = G / (np.dot(j,G@j))         # -> noramlizes for rank 1 

print(np.linalg.det(G))

D = np.linalg.det(G)

# complexity = temp_complexity((np.eye(n),s),G)

# print(complexity)


# # normalize such that JI = 1 ??

# exit()

# exit()
# G = G / (np.linalg.det(G)**(1/n)) 

# G = G / (len(s)) 

print(np.sqrt(np.linalg.det(G)))
print(G)

# exit()
# exit()
l = []
# exit()

minbad = 100000


for m1, b1 in Pmaps(search_range, s):
	# print(m1)
	badness, complexity, error = temp_measures((m1,s),G)

	# print(m1, badness, complexity)
	if badness <= 0.4: # 0.4
		minbad = badness
		# print(badness, complexity * error)
		print(m1, badness, complexity)
		# print(complexity)
		l.append((m1.copy(), badness))

l.sort(key = lambda e: e[1])
# exit()
l = l[0:20]
# l = l[-20:]

print("============================")

for t in l:
	print(t[0][0,0], t[1])

# exit()
print(len(l))
# exit()
l2 = []

r = 2

for c in itertools.combinations(l,r):
	temp = np.vstack([x[0] for x in c])
	badness, complexity, error = temp_measures((temp,s),G)
	if complexity > 0.01:
		l2.append((temp, badness, complexity))

l2.sort(key = lambda e: e[1])
l2 = l2[0:100]
for c in l2:
	# if c[1] < 3.0 and c[2] < 40:
	# if c[1] < 3.0:
	# print(c[0][0,0])
	print(c[0][:,0])
	print(hnf(c[0]))
	# print(c[1]*(2**r), c[2]*(2**r))
	print(c[1], c[2])