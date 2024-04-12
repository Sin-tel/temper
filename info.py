# utilities for the web app

import time
from typing import Optional, Any, Sequence
from markupsafe import Markup

import numpy as np
from parse import *
from lib_temper import *
from process_names import load_names, write_names

# TODO: we just stupidly do this on startup, if it starts taking to long do an actual cache
write_names()
temperament_names = load_names()


def from_commas(comma_str: str, basis: IntMat, s_expanded: SubgroupInt) -> tuple[IntMat, IntMat]:
    commas = parse_intervals(comma_str, basis, s_expanded)
    assert len(commas) > 0, "need at least one comma"
    commas = np.hstack(commas)
    commas = hnf(commas.T, remove_zeros=True).T  # fix redundant commas

    M_expanded = cokernel(commas)

    assert M_expanded[0][0] != 0, "Can't temper out the octave."

    commas_2 = solve_diophantine(basis, commas)

    M = cokernel(commas_2)

    assert np.allclose(M_expanded @ basis @ commas_2, 0), "comma not in subgroup"

    return (M, M_expanded)


def from_edos(edo_str: str, basis: IntMat, s_expanded: SubgroupInt) -> tuple[IntMat, IntMat]:
    edos = parse_edos(edo_str, get_subgroup(basis, s_expanded))

    M = hnf(np.vstack(edos), remove_zeros=True)

    # remove contorsion
    M = defactored_hnf(M)

    # find expansion with the same kernel
    M_expanded = hnf(cokernel(basis @ kernel(M)))

    return (M, M_expanded)


def search_families(s_expanded, T_expanded):
    # find temperament family names
    families: list[str] = []
    families_weak: list[str] = []

    for restrict, family_list in temperament_names.items():
        indices = subgroup_index(s_expanded, restrict)
        if indices is not None:
            r_ind = T_expanded[:, indices]
            r_ind = hnf(r_ind, remove_zeros=True)

            if r_ind.shape[0] >= r_ind.shape[1]:
                # if rank >= dim there is no point
                continue

            r_ind_weak = defactored_hnf(r_ind)

            family_index = tuple(r_ind.flatten())
            family_index_weak = tuple(r_ind_weak.flatten())
            # print(f"{restrict}, {family_index}")
            family = family_list.get(family_index)
            family_weak = family_list.get(family_index_weak)
            if family is not None:
                families.append(family)
            if family_weak is not None:
                families_weak.append(family_weak)

    families = set(families)  # type: ignore[assignment]
    families_weak = set(families_weak).difference(families)  # type: ignore[assignment]

    return families, families_weak


def subgroup_index(s: Sequence[int], l: Sequence[int]) -> Optional[list[int]]:
    res = [-1] * len(l)
    for i, k in enumerate(l):
        for j, v in enumerate(s):
            if v == k:
                res[i] = j
    for k in res:
        if k == -1:
            return None
    return res


# fmt: off
alternating_iter = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8, 9, -9, 10, -10, 11, -11, 12, -12]
# fmt: on


def info(
    T: IntMat, basis: IntMat, T_expanded: IntMat, s_expanded: SubgroupInt, options: dict[str, Any]
) -> dict[str, Any]:
    s = get_subgroup(basis, s_expanded)

    res: dict[str, Any] = {}

    # need to keep this here to get the layout correct
    # it's stupid but whatever
    res["rank"] = T.shape[0]

    res["subgroup"] = ".".join(map(str, s))

    families, families_weak = search_families(s_expanded, T_expanded)

    families_str = ""
    if len(families) > 0:
        names = [page_with_link([f"{s} family", s], s) for s in sorted(list(families))]
        families_str += ", ".join(names)
    if len(families_weak) > 0:
        if families_str != "":
            families_str += ", "
        names = [page_with_link([f"{s} family", s], s) for s in sorted(list(families_weak))]
        families_str += "(" + ", ".join(names) + ")"

    if families_str != "":
        res["families"] = families_str

    # The metric used for complexity calculations,
    # which should be the same regardless of the weights for optimization
    # (used for calculating comma and basis)
    # --
    # https://en.xen.wiki/w/Generalized_Tenney_norms_and_Tp_interval_space
    # b.T @ W^2 @ b

    G_compl_exp = np.linalg.inv(metric_wilson(s_expanded))
    G_compl = (basis.T @ G_compl_exp @ basis).astype(np.float64)

    commas = LLL(kernel(T), G_compl)

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
        vec = Markup(f'<span style="font-family: var(--mono)">[{vec}]</span>&nbsp; ')
        rat = ratio(c, s)
        comma_str.append((vec + ratio_with_link(rat)))

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
    elif T.shape[0] == 1:
        res["edo"] = edo_map_notation(T[0], s)

    # find norm on quotient space
    # WL = metric_weil_k(s_expanded, 500.0)
    # WL = T @ WL @ T.T
    # WL = np.linalg.inv(WL)
    # gens2 = LLL(np.eye(T.shape[0], dtype=np.int64), WL)
    # T = solve_diophantine(gens2, T)

    gens = preimage(T)
    gens = simplify(gens, commas, G_compl)

    # make positive
    for i in range(T.shape[0]):
        if log_interval(gens[:, i], s) < 0:
            T[i, :] = -T[i, :]
            gens[:, i] = -gens[:, i]

    if "reduce" in options:
        red = options["reduce"]
    else:
        red = "on"

    if red == "on":
        o = T[0, 0]
        genoct = np.zeros_like(gens[:, 0])
        genoct[0] = 1
        eq = log_interval(genoct, s)
        # reduce by equave
        for i in range(1, T.shape[0]):
            # replacing floor by round will find smallest
            # nb numpy floor rounds down for negative numbers, not towards zero
            red = int(np.floor(log_interval(gens[:, i], s) / eq))
            gens[0, i] -= red
            T[0, :] += o * red * T[i, :]
    elif red == "spine":
        # TODO: do the size calculations with tempered intervals instead
        o = T[0, 0]
        genoct = np.zeros_like(gens[:, 0])
        genoct[0] = 1
        eq = log_interval(genoct, s)
        # reduce first two by equave
        for i in range(1, min(2, T.shape[0])):
            red = int(np.floor(log_interval(gens[:, i], s) / eq))
            gens[0, i] -= red
            T[0, :] += o * red * T[i, :]
        if T.shape[0] > 2:
            # generally not a fifth, just the second generator
            fifth = log_interval(gens[:, 1], s)
            # reduce others
            cutoff = 0.04736875252  # Same cutoff as in modified FJS = log2(sqrt(2187/2048))
            # cutoff = 0.02781874387  # Neutral FJS = log2(sqrt(134217728/129140163))
            for i in range(2, T.shape[0]):
                interval_size = log_interval(gens[:, i], s)

                for k in alternating_iter:
                    try_next = interval_size + k * fifth
                    eq_red = -int(np.round(try_next / eq))
                    try_next += eq * eq_red
                    # print(k, eq_red, try_next)
                    if np.abs(try_next) <= cutoff:
                        gens[:, i] += k * gens[:, 1] + eq_red * gens[:, 0]
                        T[1, :] -= k * T[i, :]
                        T[0, :] -= eq_red * T[i, :]
                        break
            # we have to do this again because the resulting intervals may be complicated
            gens = simplify(gens, commas, G_compl)
            # make positive
            for i in range(T.shape[0]):
                if log_interval(gens[:, i], s) < 0:
                    T[i, :] = -T[i, :]
                    gens[:, i] = -gens[:, i]

    res["mapping"] = format_matrix(T)

    gens_print = [ratio_with_link(ratio(g, s)) for g in gens.T]

    res["preimage"] = gens_print

    if "weights" in options:
        weight = options["weights"]
    else:
        weight = "weil"

    # legacy
    if options["tenney"]:
        weight = "tenney"

    j2 = log_subgroup(s)

    # get the equave
    equave = factors(s[0], s_expanded)

    te_tun, _ = lstsq((T_expanded, s_expanded), weight)
    cte_tun, _ = cte((T_expanded, s_expanded), weight, V=equave)

    showtarget = False
    if "target" in options:
        targets = parse_intervals(options["target"], basis, s_expanded)
        showtarget = len(targets) > 0

    if showtarget:
        targets = np.hstack(targets)
        n_targets = targets.shape[1]

        # use least squares if theres more targets than degrees of freedom
        # if equal, this just solves the system
        if n_targets >= T_expanded.shape[0]:
            target_tun, _ = lstsq((T_expanded, s_expanded), weight="unweighted", V=targets)

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

    W_tenney = basis.T @ np.linalg.inv(metric_tenney(s_expanded)) @ basis
    W_tenney = W_tenney.astype(np.float64)
    W_tenney_inv = np.linalg.inv(W_tenney)

    badness = temp_badness((T, s), W=W_tenney_inv)
    res["badness"] = f"{badness:.3f}"
    # complexity = temp_complexity((T_expanded, s_expanded))
    # res["complexity"] = f"{complexity:.2f}"

    return res
