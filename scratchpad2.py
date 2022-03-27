from util import *
import math
from search_1_53__2_3_5 import templist

def makemap(t,r):
	return np.reshape(np.array(t),(r,-1))

S = [2,3,5]

nps = np.array(S)

def wilson(c):
	return np.sum(np.abs(c.T) * nps)

list2 = []

j = np.log2(nps)
W = np.diag(1. / j)

# for i in range(2000):
for i in range(len(templist)):
	# print("===========================================")
	if i%100 == 0:
		print(i/len(templist))

	M = makemap(templist[i],2)
	c = kernel(M)

	c = make_positive(c,S)

	r, d = M.shape

	# print(r,d)

	# compl = wilson(c)*(0.015+log_interval(c, S))

	sol, e = lstsq((M,S), weight="tenney")

	# BREED
	# err = np.sqrt(np.average((e @ W)**2)) # absolute error
	# compl  = np.sqrt (np.linalg.det ((M @ W) @ (M @ W).T / M.shape[1])) # normalize by d^r

	# SMITH
	err = np.sqrt(np.sum((e @ W)**2) * (r + 1) / (d-r) )
	compl = np.sqrt (np.linalg.det ((M @ W) @ (M @ W).T ) / math.comb(d,r))

	# print(err, compl)

	badness = err*(compl**(d/(d-r)))  # logflat

	# if compl < 5:
	list2.append([M,c,badness,err,compl])  
	# list2.append([M,c,compl])

list2.sort(key = lambda l:l[2])

for k in list2:
	print("===========")
	M = k[0]
	print(M)
	# print(M @ W @ W @ M.T)
	print(k[1].T)
	print(ratio(k[1],S))
	print(k[2])
	# print(k[3]*1200)
	# print(k[4])


