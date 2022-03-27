from util import *

# from itertools import combinations

# plist = list(combinations([1,2,3,4,5,6,7,8,9,10], 3))

# plist.sort(key = sum)
# print(plist)

S = p_limit(11)


logs = np.log2(np.array(S))

search_range = (0.5,3110.5)

# biglist = []

# schisma = factors("32805/32768", S)
syntonic = factors("81/80", S)
starling = factors("126/125", S)
c3 = factors("99/98", S)
# c4 = factors("66/65", S)
# c5 = factors("51/50", S)
# c6 = factors("57/56", S)

# c = np.hstack([syntonic, starling, c3, c4, c5, c6])

# c = np.hstack([schisma])
c = np.hstack([syntonic,starling,c3])
# print(c)

# print(hnf(np.array([[2,3,4],[3,5,8]])))

e_max = 0.08

for m1,b1 in Pmaps(search_range, S):
	if np.all(m1 @ c == 0):

		tun = np.sum((m1/logs)**2)/np.sum(m1/logs)  # optimal TE tuning for edo

		e = np.sqrt(np.average((m1/logs - tun)**2))  # relative error

		# e2 = np.sqrt(np.average((err/logs)**2))*tun  # same



		# print(tun,tun2,tun3)

		# e = e * tun ** (1.0 / (len(S) - 1)) # logflat

		if e < e_max:
			# e_max = e

			print(tun, e) # RMSe

			# print(m1, tun, np.round(e*1000)/1000)
	# if m1 @ syntonic == 0:
		# print(m1)