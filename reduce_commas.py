from lib_temper import *
import numpy as np

# from hsnf import column_style_hermite_normal_form, row_style_hermite_normal_form, smith_normal_form
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

subgroup = p_limit(53)
# subgroup = p_limit(19)
# subgroup = p_limit(7)

j = log_subgroup(subgroup)
d = j.shape[0]

# edo = patent_map(12, subgroup)
edo = patent_map(20567, subgroup)
print(edo)


def kernel(M):
    assert M.ndim == 2
    r, d = M.shape

    M = np.vstack([M, np.eye(d, dtype=np.int64)])
    print(M)
    print(Matrix(M.T))
    H = hermite_normal_form(Matrix(M.T))
    print(np.array(H))

    K = np.array(H).T
    print(K)
    K = K[r::, r::]


#     return K


commas = kernel(edo)

print(edo @ commas)

exit()

W = np.linalg.inv(metric_wilson(subgroup))
# W = np.linalg.inv(metric_weil(subgroup))
# W = np.linalg.inv(metric_tenney(subgroup))
# W = np.eye(d)

commas = LLL(commas, W)
print(commas.T)
for i in range(commas.shape[1]):
    print(ratio(commas[:, i], subgroup))
