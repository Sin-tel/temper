# This code is just copy pasted from x31eq
# http://x31eq.com/temper.html
# thanks to Graham Breed!

# Uses lists of python ints for everything
# Which might be slower than numpy, but they are guaranteed to not overflow
# In practice, this seems to be pretty fast

import numpy as np


def hnf_bigint(m):
    h = lattice_reduction(m.tolist())
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
