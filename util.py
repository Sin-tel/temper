# utilities for the web app

import numpy as np
import re
from lib_temper import *


def cents(x, prec=3):
    return "{1:.{0}f}".format(prec, 1200 * x)


def parse_subgroup(s):
    s = [Fraction(i) for i in re.split("[\\.,; ]+", s)]

    if len(s) == 1:
        expanded = p_limit(s[0])
        return np.eye(len(expanded), dtype=np.int64), expanded
    else:
        s_basis, expanded = get_subgroup_basis(s)
        # s = get_subgroup(s_basis, expanded)
        return s_basis, expanded


def parse_edos(s):
    s = [int(i) for i in re.split("[\\.,; &]+", s)]
    return s


ratioPattern = "(\d+)[/:](\d+)"
vectorPattern = "[\[(<]\s*(-?\d+(?:[,\s]+-?\d+)*)\s*[\])>]"


def parse_commas(c, s):
    commas = []
    for n, d in re.findall(ratioPattern, c):
        print(n, d, flush=True)
        commas.append(factors((int(n), int(d)), s))
    for v in re.findall(vectorPattern, c):
        l = len(s)
        res = np.zeros((l, 1), dtype=np.int64)
        v = np.array(list(map(int, v.replace(",", " ").split())))
        res[: v.shape[0], 0] = v[:l]
        commas.append(res)

    return commas


def format_matrix(matrix):
    """Format a matrix using LaTeX syntax"""
    body_lines = [" & ".join(map(str, row)) for row in matrix]
    body = "\\\\\n".join(body_lines)
    body = "\\begin{bmatrix}" + body + "\\end{bmatrix}"
    return "\\( " + body + " \\)"


def from_commas(args):
    basis, s_expanded = parse_subgroup(args["subgroup"])

    commas = parse_commas(args["commas"], s_expanded)

    commas = np.hstack(commas)

    commas = hnf(commas.T, remove_zeros=True).T  # fix redundant commas

    M_expanded = cokernel(commas)

    assert M_expanded[0][0] != 0, "Can't temper out the octave."

    commas_2 = solve_diophantine(basis, commas)

    M = cokernel(commas_2)

    assert np.allclose(M_expanded @ basis @ commas_2, 0), "comma not in basis"

    return (M, basis, M_expanded, s_expanded)


def from_edos(args):
    basis, s_expanded = parse_subgroup(args["subgroup"])
    edo_list = parse_edos(args["edos"])

    edos = []
    for e in edo_list:
        edos.append(patent_map(e, get_subgroup(basis, s_expanded)))

    M = hnf(np.vstack(edos), remove_zeros=True)

    # remove contorsion
    M = defactored_hnf(M)

    # find expansion with the same kernel
    M_expanded = hnf(cokernel(basis @ kernel(M)))

    return (M, basis, M_expanded, s_expanded)


# given a map for some edo, return the string representation in the new format
# (to replace warts notation)
def edo_map_notation(this_map, subgroup):
    this_edo = this_map[0]
    j = log_subgroup(subgroup)

    patent_map = np.round(this_edo * j).astype(np.int64)
    diff = this_map - patent_map
    adjustments = []
    for i, p in enumerate(diff):
        if p != 0:
            sign = ""
            if p < 0:
                sign = "-" * abs(p)
            elif p > 0:
                sign = "+" * p
            adjustments.append(sign + str(subgroup[i]))

    mstr = str(this_edo)
    if len(adjustments) > 0:
        mstr += "[" + ", ".join(adjustments) + "]"
    return mstr


def info(temp, options):
    T = temp[0]
    basis = temp[1]
    T_expanded = temp[2]
    s_expanded = temp[3]

    s = get_subgroup(basis, s_expanded)

    res = dict()

    res["rank"] = T.shape[0]
    res["dim"] = T.shape[1]

    res["subgroup"] = ".".join(map(str, s))

    s_w = []
    for fr in s:
        fr = fr.as_integer_ratio()
        s_w.append(fr[0] * fr[1])
    W_wilson = np.diag(s_w).astype(np.double)

    G_wilson = W_wilson @ W_wilson.T

    commas = LLL(kernel(T), G_wilson)

    for i in range(len(commas.T)):
        commas[:, i] = make_positive(commas[:, i], s)

    # G_wilson_dual = np.linalg.inv(G_wilson)

    c_length = []
    for c_column in commas:
        c_length.append(max([len(str(x)) + 1 for x in list(c_column)]))

    comma_str = []
    for c in commas.T:
        ## format spaces
        vec = "".join(["{1: {0}d}".format(c_length[i], c[i]) for i in range(len(c))])
        vec = vec.replace(" ", "&nbsp;")
        vec = (
            '<span style="font-family: var(--mono)">['
            + vec
            + "]<sup>T</sup></span>&nbsp; "
        )
        comma_str.append(vec + str(ratio(c, s)))

    res["comma basis"] = "<br>".join(comma_str)

    edolist = find_edos(T, s)
    # edolist = find_edos_patent(T,s)

    if edolist is not None and len(edolist) >= 1:
        maps_join = find_join(T, s, edolist)

        show_list = []
        for m in edolist:
            mstr = edo_map_notation(m[0][0], s)
            show_list.append(mstr)

        res["edos"] = ", ".join(map(str, show_list))

        if maps_join is not None:
            res["edo join"] = " & ".join(
                map(lambda x: edo_map_notation(x, s), maps_join)
            )

    # T = LLL(T.T, G_wilson_dual).T

    gens = preimage(T)

    if options["reduce"]:
        # eq = log_interval(gens[0], s)
        o = T[0, 0]
        genoct = np.zeros_like(gens[:, 0])
        genoct[0] = 1
        # reduce by octave
        for i in range(1, T.shape[0]):
            # make positive first
            if log_interval(gens[:, i], s) < 0:
                T[i, :] = -T[i, :]
                gens[:, i] = -gens[:, i]

            red = int(np.floor(log_interval(gens[:, i], s)))
            gens[:, i] -= red * genoct
            T[0, :] += o * red * T[i, :]

    else:
        # make positive
        for i in range(T.shape[0]):
            if log_interval(gens[:, i], s) < 0:
                T[i, :] = -T[i, :]
                gens[:, i] = -gens[:, i]

    gens = simplify(gens, commas, G_wilson)

    res["mapping"] = format_matrix(T)

    gens_print = [ratio(g, s) for g in gens.T]

    res["preimage"] = list(map(str, gens_print))

    weight = "unweighted"
    if options["tenney"]:
        weight = "tenney"

    j2 = log_subgroup(s)

    te_tun, te_err = lstsq((T_expanded, s_expanded), weight)
    cte_tun, cte_err = cte((T_expanded, s_expanded), weight)

    showtarget = False
    if "target" in options:
        targets = parse_commas(options["target"], s_expanded)
        showtarget = len(targets) > 0

    if showtarget:
        targets = np.hstack(targets)
        n_targets = targets.shape[1]

        # use least squares if theres more targets than constraints
        # (or equal, which will just solve the system)
        # for soem reason it only works unweighted, need to investigate
        if n_targets >= T_expanded.shape[0]:
            target_tun, target_err = lstsq(
                (T_expanded, s_expanded), weight="unweighted", V=targets
            )

        # otherwise, use constraints with reguular weights
        else:
            target_tun, target_err = cte((T_expanded, s_expanded), weight, V=targets)

        target_tun2 = (target_tun.T @ T_expanded @ basis) @ gens

        target_err2 = target_tun2 @ T - j2

    te_tun2 = (te_tun.T @ T_expanded @ basis) @ gens
    cte_tun2 = (cte_tun.T @ T_expanded @ basis) @ gens

    print(te_tun.T @ T_expanded)
    print((te_tun.T @ T_expanded @ basis) @ gens @ T)

    te_err2 = te_tun2 @ T - j2
    cte_err2 = cte_tun2 @ T - j2

    res["TE tuning"] = list(map(cents, te_tun2.flatten()))
    res["CTE tuning"] = list(map(cents, cte_tun2.flatten()))
    if showtarget:
        target_str = []
        for c in targets.T:
            target_str.append(str(ratio(c, s_expanded)))
        res["target tuning (" + ", ".join(target_str) + ")"] = list(
            map(cents, target_tun2.flatten())
        )

    res["TE errors"] = ", ".join(map(cents, te_err2.flatten()))
    res["CTE errors"] = ", ".join(map(cents, cte_err2.flatten()))
    if showtarget:
        res["target errors"] = ", ".join(map(cents, target_err2.flatten()))

    # print((cte_tun.T @ T_expanded @ basis))
    # res["CTE map"] = " ".join(map(cents, (cte_tun.T @ T_expanded @ basis)))

    # badness, complexity, error = temp_measures((T_expanded, s_expanded))
    # res["badness"] = '{:.3e}'.format(badness)
    # res["complexity"] = '{1:.{0}f}'.format(3, complexity)
    # res["error"] = '{1:.{0}f}'.format(3, 1200*error)

    return res


#########################################
if __name__ == "__main__":
    args = dict()
    # args["subgroup"] = "2.3.5.7"
    args["subgroup"] = "2.5.9/7"
    # args["edos"] = "19,22"
    args["commas"] = "225/224"

    options = dict()

    options["tenney"] = True
    options["reduce"] = True

    # temp = from_edos(args)
    temp = from_commas(args)
    html_info = info(temp, options)
