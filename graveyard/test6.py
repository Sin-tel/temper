from util import *
import numpy as np
import math


args = dict()

args["commas"] = "3025/3024 6000/5929"
args["subgroup"] = "2.5/3.7/3.11/3"

# args["commas"] = "81/80"
# args["subgroup"] = "2.3/2.5"


M, basis, M_expanded, s_expanded = from_commas(args)

s = get_subgroup(basis, s_expanded)

j = log_subgroup(s_expanded)

W_e = np.diag(1 / j)


W = W_e @ basis

print(W)

print(basis)
print(W_e @ W_e.T)
print(W @ W.T)
# print(basis.T @ W_e)

gens = preimage(M)

j = np.atleast_2d(j)

G_e = W_e.T @ W_e
# sol = np.linalg.lstsq((M_expanded @ W_e).T, (j @ W_e).T, rcond=None)[0]
sol = np.linalg.inv(M_expanded @ G_e @ M_expanded.T) @ M_expanded @ G_e @ j.T

tun = sol.T @ M_expanded
te_err = (tun - j).flatten()

te_tun = sol

print("================")
print(1200 * te_tun)


te_tun2 = (te_tun.T @ M_expanded @ basis) @ gens

print(1200 * te_tun2)

# print(M_expanded @ W_e)

# print(M)
# print(W.T @ W)


j2 = log_subgroup(s)
j2 = np.atleast_2d(j2)
# te_tun3 = np.linalg.lstsq((M @ W.T).T, (j2 @ W.T).T, rcond=None)[0]

G = W.T @ W

G_inv = np.linalg.inv(G_e)
G2 = np.linalg.inv(basis.T @ G_inv @ basis)

# print(G)
te_tun3 = np.linalg.inv(M @ G2 @ M.T) @ M @ G2 @ j2.T

# print(1200*te_tun3.T)

print("===============")

print(G2)

print(basis.T @ G_inv @ basis)


# print(j @ basis)
# print(j2)

# l = np.array([[-1,-1,2,0]]).T

# l_e = basis @ l

# print((W_e @ l_e).T @ (W_e @ l_e))
# print((W @ l).T @ (W @ l))


# print(M_expanded)

# print(M_expanded @ basis)

# print(basis.T @ G_e @ basis)

# print(G_e)
# print(G)
