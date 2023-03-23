from util import *

# S = p_limit(7)
# c = factors("225/224",S)
# print(c)
# M = cokernel(np.hstack([c,factors("2/1",S)]))
# print(M)
# Mr = LLL(M.T).T

# print(Mr)
# pre = (preimage(Mr))
# for v in pre:
# 	print(ratio(v,S))


S = p_limit(5)

search_range = (0.5, 53.5)

biglist = []

for m1, b1 in Pmaps(search_range, S):
    i1 = (b1[0] + b1[1]) / 2.0
    print(i1)
    for m2, b2 in Pmaps(search_range, S):
        if b1 > b2:
            i2 = (b2[0] + b2[1]) / 2.0
            M = np.vstack([m1, m2])
            if factor_order(M) == 1:
                M = hnf(M)
                # print("=====")
                # print(i1,i2)
                # print(M)
                biglist.append(tuple(M.flatten()))

biglist = sorted(list(set(biglist)))

for t in biglist:
    print(t)

print(len(biglist))
