import numpy as np
from lib_temper import *


def hnf_bigint(m):
    h = lattice_reduction(m)
    return np.array(h, dtype=np.int64)


def lattice_reduction(ets):
    """Reduced echelon form without removing common factors
    and allowing positive numbers to appear in place of zeros
    instead of introducing torsion.
    (These positive numbers should be the smallest possible.)
    """
    echelon = echelon_form(ets)
    n_columns = len(echelon[0])
    column = 1
    for row in range(1, len(echelon)):
        nrow = echelon[row]
        try:
            while nrow[column] == 0:
                column += 1
        except IndexError:
            # lower rank
            break
        n = nrow[column]
        for srow in echelon[:row]:
            m = srow[column] // n
            for j in range(n_columns):
                srow[j] -= m * nrow[j]
    return echelon


def echelon_form(ets, col=0):
    working = [normalize_sign(et) for et in map(list, ets)]
    if not working:
        return []
    ncols = len(working[0])
    if col == ncols:
        return working
    assert col < ncols

    reduced = []
    while True:
        reduced += [row for row in working if row[col] == 0]
        working = [row for row in working if row[col] != 0]

        if len(working) < 2:
            return working + echelon_form(reduced, col + 1)

        working.sort()
        pivot = working[0]
        for row in working[1:]:
            n = row[col] // pivot[col]
            for i in range(col, ncols):
                row[i] -= pivot[i] * n


def normalize_sign(mapping):
    nonzeros = [m for m in mapping if m]
    if nonzeros == [] or nonzeros[0] > 0:
        return mapping
    return [-m for m in mapping]


subgroup = [2, 3, 5, 7]


M = np.vstack(
    [
        patent_map(12, subgroup),
        patent_map(22, subgroup),
        patent_map(31, subgroup),
        patent_map(53, subgroup),
    ]
)

print(M)
print(M.tolist())
print(hnf_bigint(M))
print(hnf(M))
