# utilities for the web app

import re
import sys
import time

import numpy as np
from lib_temper import *


def cents(x, prec=3):
    return "{1:.{0}f}".format(prec, 1200 * x)


def parse_subgroup(s):
    s = [Fraction(i) for i in re.split(r"[\\.,; ]+", s)]

    if len(s) == 1:
        expanded = p_limit(s[0])
        return np.eye(len(expanded), dtype=np.int64), expanded
    else:
        s_basis, expanded = get_subgroup_basis(s)
        return s_basis, expanded


wart_map = {
    "a": 2,
    "b": 3,
    "c": 5,
    "d": 7,
    "e": 11,
    "f": 13,
    "g": 17,
    "h": 19,
    "i": 23,
    "j": 29,
    "k": 31,
    "l": 37,
    "m": 41,
    "n": 43,
    "o": 47,
}


def parse_edos(s, subgroup):
    # split string by separators (,.;& ) but ignore when in square brackets
    # "12, 17[+5, 8], 22"
    # => ["12", "17[+5, 8]", "22"]
    s = re.split(r"[\\.,; &]+(?!(?:[^,\[\]]+,)*[^,\[\]]+])", s.lower().strip())

    edos = []

    log_s = log_subgroup(subgroup)
    log_s = log_s / log_s[0]  # fix equave

    for e in s:
        res = re.split(r"(\d+)", e)[1:]
        if res == []:
            continue
        edo_num = int(res[0])
        p_map = patent_map(edo_num, subgroup)
        if res[1] == "":
            # if the input was simply an integer, then add the patent map
            edos.append(p_map)
        else:
            adjust = re.findall(r"\[.*\]", e)
            if len(adjust) > 0:
                # here we try to parse formats like '17[+5]'
                adjust = adjust[0][1:-1]
                adjust = re.split(r"[\\.,; &]+", adjust)
                for l in adjust:
                    p = re.split(r"([+-]+)", l)[1:]
                    adj_str = p[0]
                    ratio = re.split(r"[/:]", p[1])
                    if len(ratio) == 1:
                        ratio = Fraction(int(ratio[0]))
                    else:
                        ratio = Fraction(int(ratio[0]), int(ratio[1]))

                    index = None
                    for i, pr in enumerate(subgroup):
                        if ratio == pr:
                            index = i
                            break
                    if index is None:
                        raise AssertionError("Adjustment not in subgroup")

                    for c in adj_str:
                        if c == "+":
                            p_map[0, index] = p_map[0, index] + 1
                        elif c == "-":
                            p_map[0, index] = p_map[0, index] - 1

                edos.append(p_map)
            else:
                # here we try to parse formats like '17c' (aka wart notation)
                warts = re.finditer(r"(\D)\1*", res[1])
                for g in warts:
                    w = g.group()
                    w_count = len(w)
                    # p means patent so we don't have to  adjust anything
                    if w[0] != "p":
                        w_prime = wart_map[w[0]]

                        index = None
                        for i, p in enumerate(subgroup):
                            assert p.denominator == 1, "Warts can't be used in rational subgroups"
                            if w_prime == p.numerator:
                                index = i
                                break
                        if index is None:
                            raise AssertionError("Wart not in subgroup")

                        # the actual adjustment count is only half the number of letters
                        count = (w_count + 1) // 2

                        float_prime = edo_num * log_s[index]
                        sign = 1
                        # check if the first one has to round up or down
                        if float_prime - np.floor(float_prime + 0.5) <= 0:
                            sign *= -1
                        # alternate added or subtracting
                        if w_count % 2 == 0:
                            sign *= -1

                        p_map[0, index] = p_map[0, index] + sign * count

                edos.append(p_map)

    return edos


ratio_pattern = r"(\d+)[/:](\d+)"
vector_pattern = r"[\[(<]\s*(-?\d+(?:[,\s]+-?\d+)*)\s*[\])>]"


def parse_intervals(c, s):
    commas = []
    for n, d in re.findall(ratio_pattern, c):
        commas.append(factors((int(n), int(d)), s))
    for v in re.findall(vector_pattern, c):
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

    commas = parse_intervals(args["commas"], s_expanded)

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
    edos = parse_edos(args["edos"], get_subgroup(basis, s_expanded))

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

    # divide by equave to fix non-octave temps
    j = log_subgroup(subgroup) / np.log2(float(subgroup[0]))

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


def name_tuning_weight(weight):
    if weight == "unweighted":
        return "E"
    elif weight == "tenney":
        return "TE"
    elif weight == "weil":
        return "WE"
    else:
        raise ValueError("unknown weight parameter")


def info(temp, options):
    T = temp[0]
    basis = temp[1]
    T_expanded = temp[2]
    s_expanded = temp[3]

    s = get_subgroup(basis, s_expanded)

    res = {}

    # need to keep this here to get the layout correct
    # it's stupid but whatever
    res["rank"] = T.shape[0]

    res["subgroup"] = ".".join(map(str, s))

    # https://en.xen.wiki/w/Generalized_Tenney_norms_and_Tp_interval_space
    # b.T @ W^2 @ b
    G_wilson_exp = np.linalg.inv(metric_wilson(s_expanded))
    G_wilson = basis.T @ G_wilson_exp @ basis

    commas = LLL(kernel(T), G_wilson)

    for i in range(len(commas.T)):
        commas[:, i] = make_positive(commas[:, i], s)

    c_length = []
    for c_column in commas:
        c_length.append(max(len(str(x)) + 1 for x in list(c_column)))

    comma_str = []
    for c in commas.T:
        ## format spaces
        vec = "".join(["{1: {0}d}".format(c_length[i], c[i]) for i in range(len(c))])
        vec = vec.replace(" ", "&nbsp;")
        vec = '<span style="font-family: var(--mono)">[' + vec + "]<sup>T</sup></span>&nbsp; "
        comma_str.append(vec + str(ratio(c, s)))

    res["comma basis"] = "<br>".join(comma_str)

    t_start = time.time()

    edolist = find_edos(T, s)

    print("find_edos() took: ", time.time() - t_start)

    if edolist is not None and len(edolist) >= 1:
        maps_join = find_join(T, s, edolist)

        show_list = []
        for m in edolist:
            mstr = edo_map_notation(m[0][0], s)
            show_list.append(mstr)

        res["edos"] = ", ".join(map(str, show_list))

        if maps_join is not None:
            res["edo join"] = " & ".join(map(lambda x: edo_map_notation(x, s), maps_join))

    # find norm on quotient space
    # WL = metric_weil_k(s_expanded, 10.0)
    # WL = T @ WL @ T.T
    # WL = np.linalg.inv(WL)

    # gens2 = LLL(np.eye(T.shape[0], dtype=np.int64), WL)
    # T = solve_diophantine(gens2, T)

    gens = preimage(T)
    gens = simplify(gens, commas, G_wilson)

    if options["reduce"]:
        o = T[0, 0]
        genoct = np.zeros_like(gens[:, 0])
        genoct[0] = 1
        # eq = log_interval(gens[:, 0], s)
        eq = log_interval(genoct, s)
        # reduce by equave
        for i in range(1, T.shape[0]):
            # replacing floor by round will find smallest
            # nb numpy floor rounds down for negative numbers, not towards zero
            red = int(np.round(log_interval(gens[:, i], s) / eq))
            gens[0, i] -= red
            T[0, :] += o * red * T[i, :]

    # make positive
    for i in range(T.shape[0]):
        if log_interval(gens[:, i], s) < 0:
            T[i, :] = -T[i, :]
            gens[:, i] = -gens[:, i]

    res["mapping"] = format_matrix(T)

    gens_print = [ratio(g, s) for g in gens.T]

    res["preimage"] = list(map(str, gens_print))

    if "weights" in options:
        weight = options["weights"]

    # legacy
    if options["tenney"]:
        weight = "tenney"

    j2 = log_subgroup(s)

    # get the equave
    equave = factors(s[0], s_expanded)

    te_tun, te_err = lstsq((T_expanded, s_expanded), weight)
    cte_tun, cte_err = cte((T_expanded, s_expanded), weight, V=equave)

    showtarget = False
    if "target" in options:
        targets = parse_intervals(options["target"], s_expanded)
        showtarget = len(targets) > 0

    if showtarget:
        targets = np.hstack(targets)
        n_targets = targets.shape[1]

        # use least squares if theres more targets than constraints
        # (or equal, which will just solve the system)
        # for some reason it only works unweighted, need to investigate
        if n_targets >= T_expanded.shape[0]:
            target_tun, target_err = lstsq((T_expanded, s_expanded), weight="unweighted", V=targets)

        # otherwise, use constraints with regular weights
        else:
            target_tun, target_err = cte((T_expanded, s_expanded), weight, V=targets)

        target_tun2 = (target_tun.T @ T_expanded @ basis) @ gens

        target_err2 = target_tun2 @ T - j2

    te_tun2 = (te_tun.T @ T_expanded @ basis) @ gens
    cte_tun2 = (cte_tun.T @ T_expanded @ basis) @ gens

    te_err2 = te_tun2 @ T - j2
    cte_err2 = cte_tun2 @ T - j2

    weight_name = name_tuning_weight(weight)

    res[f"{weight_name} tuning"] = list(map(cents, te_tun2.flatten()))
    res[f"C{weight_name} tuning"] = list(map(cents, cte_tun2.flatten()))
    if showtarget:
        target_str = []
        for c in targets.T:
            target_str.append(str(ratio(c, s_expanded)))
        res["target tuning (" + ", ".join(target_str) + ")"] = list(
            map(cents, target_tun2.flatten())
        )

    res[f"{weight_name} errors"] = ", ".join(map(cents, te_err2.flatten()))
    res[f"C{weight_name} errors"] = ", ".join(map(cents, cte_err2.flatten()))
    if showtarget:
        res["target errors"] = ", ".join(map(cents, target_err2.flatten()))

    badness = temp_badness((T_expanded, s_expanded))
    res["badness"] = "{:.3f}".format(badness)
    # res["complexity"] = "{1:.{0}f}".format(3, complexity)
    # res["error"] = "{1:.{0}f}".format(3, 1200 * error)

    return res
