from util import *
import math

np.random.seed(2022)

## use log error / complexity

def temp_error(temp):
	M, S = temp
	r, d = M.shape

	sol, e = lstsq(temp, weight="tenney")


	# BREED
	err = np.sqrt(np.average((e @ W)**2)) 

	# SMITH
	# err = np.sqrt(np.sum((e @ W)**2) * (r + 1) / (d-r) ) 

	return err

def temp_complexity(temp):
	M, S = temp
	r, d = M.shape

	# BREED
	compl  = np.sqrt (np.linalg.det ((M @ W) @ (M @ W).T / d))

	# SMITH
	# compl = np.sqrt (np.linalg.det ((M @ W) @ (M @ W).T ) / math.comb(d,r))

	# print((M) @ (M).T)
	# print(np.sqrt(np.linalg.det(M @ M.T))**(1/r))

	# compl = np.sqrt (np.linalg.det ((M @ W) @ (M @ W).T ) / math.comb(d,r)) ** (1/r)
	# compl = (r) * np.sqrt (np.linalg.det ((M @ W) @ (M @ W).T / M.shape[1] )) ** (1/r)
	return compl

def temp_badness(temp):
	M, S = temp
	r, d = M.shape
	# logflat
	badness = temp_error(temp)*(temp_complexity(temp)**(d/(d-r)))  

	## (1/(d-r) + 1/r) (1/r)

	return badness



def tprint(m,S):
	print("=====================================")
	print(m)
	print(1200*temp_error((m,S)))
	print(temp_complexity((m,S)))
	print(temp_badness((m,S)))

iters = 1000



for d in range(10,20):
	# print(d)
	# print(primes[d])
	S = primes[0:d]
	# print(S)

	j = np.log2(np.array(S))
	W = np.diag(1. / j)

	# for r in range(2,3):
	# for r in range(2,d):
	for r in range(9,10):
		print((d,r))

		# print(math.comb(d,r))

		etot = 0
		ctot = 0
		btot = 0

		for k in range(iters):
			l = []
			for i in range(r):
				m = np.random.rand()*295+5
				l.append(patent_map(m,S))
			l = np.vstack(l)
			# print(hnf(l))

			if np.linalg.det(l @ l.T) >= 0.5:

				e = 1200*temp_error((l,S))
				c = temp_complexity((l,S))
				b = e*(c**(d/(d-r)))

				# e = np.log(temp_error((l,S)))
				# c = np.log(temp_complexity((l,S)))
				# b2 = e + (c*(d/(d-r)))

				# print(np.log(b) - b2)
				# b = e*c

				etot += e
				ctot += c
				btot += b

		print((etot/iters))
		print(ctot/iters)
		print(btot/iters)




