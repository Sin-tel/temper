from util import *
import math

np.random.seed(2022)

## use log error / complexity


def tprint(m,S):
	print("=====================================")
	print(m)
	print(1200*temp_error((m,S)))
	print(temp_complexity((m,S)))
	print(temp_badness((m,S)))

iters = 100



for d in range(2,10):
	# print(d)
	# print(primes[d])
	S = primes[0:d]
	# print(S)

	j = np.log2(np.array(S))
	W = np.diag(1. / j)

	# for r in range(2,3):
	# for r in range(2,d):
	for r in range(1,d):
		print((d,r))

		# print(math.comb(d,r))

		etot = 0
		ctot = 0
		btot = 0

		for k in range(iters):
			l = []
			for i in range(r):
				m = np.random.rand()*53+5
				l.append(patent_map(m,S))
			l = np.vstack(l)
			# print(hnf(l))

			if np.linalg.det(l @ l.T) >= 0.5:

				e = temp_error((l,S))
				c = temp_complexity((l,S))
				b = e*(c**(d/(d-r)))

				etot += e
				ctot += c
				btot += b

		# print((etot/iters))
		# print(ctot/iters)
		print(btot/iters)




