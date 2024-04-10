from typing import Any
import time
import numpy as np
import flask

from info import from_commas, search_families
from lib_temper import *
from parse import *

np.set_printoptions(suppress=True)


def temperament_search(args: dict[str, Any]) -> dict[str, Any]:
    basis, s_expanded = parse_subgroup(args["subgroup"])

    s = get_subgroup(basis, s_expanded)

    # check if we have a normal prime-basis subgroup
    s_non_prime = True
    if (basis.shape[0] == basis.shape[1]) and (basis == np.eye(basis.shape[0])).all():
        s_non_prime = False

    badness_type = "cangwu"
    if "badness" in args:
        if args["badness"] == "dirichlet":
            badness_type = "dirichlet"

    t_map = np.eye(len(s), dtype=np.int64)
    if "commas" in args:
        t_map, _ = from_commas(args["commas"], basis, s_expanded)

    res: dict[str, Any] = {}
    # search down

    rank = t_map.shape[0]
    dim = t_map.shape[1]

    W_tenney = basis.T @ np.linalg.inv(metric_tenney(s_expanded)) @ basis
    W_tenney = W_tenney.astype(np.float64)
    W_tenney_inv = np.linalg.inv(W_tenney)

    W_wilson = basis.T @ np.linalg.inv(metric_wilson(s_expanded)) @ basis
    W_wilson = W_wilson.astype(np.float64)

    W_weil = basis.T @ np.linalg.inv(metric_weil_k(s_expanded, 800.0)) @ basis
    W_weil = W_weil.astype(np.float64)
    W_weil_inv = np.linalg.inv(W_weil)

    # normalize W so it has det(W) = 1
    W_det = np.linalg.det(W_weil_inv)
    assert W_det > 0
    W_weil_inv = W_weil_inv / np.power(W_det, 1 / W_weil_inv.shape[0])

    log_s = np.log2(np.array(s).astype(np.float64))

    f_init = 32 * rank * np.max(log_s)
    f = 1.4 * (1.2 ** (dim - 1))
    # f = 2.0

    print(f"f init = {f_init}")
    print(f"f = {f}")

    eq = log_s[0]
    log_s = log_s[1:]
    log_s = np.atleast_2d(log_s)

    B = np.block([[np.eye(dim)], [f_init * log_s.T, -f_init * eq * np.eye(dim - 1)]])
    # B = np.block(
    #     [[np.eye(dim)], [f_init * np.ones((dim - 1, 1)), -f_init * eq * np.diagflat(1 / log_s)]]
    # )
    # print(B)
    B = B @ t_map.astype(np.float64).T
    # print(B)

    W_LLL = np.block(
        [[W_tenney_inv, np.zeros((dim, dim - 1))], [np.zeros((dim - 1, dim)), np.eye(dim - 1)]]
    )

    # print(B)
    # print(W_LLL)

    by_rank: dict[int, list[tuple[str, Any, float]]] = {}
    for r_s in range(rank - 1, 0, -1):
        by_rank[r_s] = []

    t_start = time.time()
    total_count = 0

    checked: set[tuple[int, ...]] = set()

    f_prev = None

    for i in range(8):
        if total_count > 100:
            break
        B = LLL(B, W=W_LLL, delta=0.9)
        # B = LLL(B, delta=0.9)
        # B = LLL(B, delta=0.9)

        B[dim:] *= f

        found = np.abs(B[0:dim].T)

        if np.any(found[:, 0] > 1000):
            break

        assert np.allclose(found - np.round(found), 0)
        found = np.round(found).astype(np.int64)

        if f_prev and f_prev == tuple(found[:, 0]):
            break
        f_prev = tuple(found[:, 0])
        print(f_prev)

        edo_list = []
        for k in found:
            if 2 < k[0] <= 1000:
                t_edo = k[None, :]
                label = edo_map_notation(k, s)
                if badness_type == "dirichlet":
                    b = temp_badness((k[None, :], s), W_tenney_inv)
                else:
                    b = height(k[None, :], W_weil_inv)

                r_edo = (label, k, b)

                t_tup = tuple(t_edo.flatten())
                if t_tup not in checked:
                    by_rank[1].append(r_edo)
                checked.add(t_tup)

                edo_list.append(r_edo)
        edo_list = sorted(edo_list, key=lambda v: v[2])

        for r_s in range(2, rank):
            count = 0
            # for combo in itertools.combinations(edo_list, r_s):
            for idx in comboBySum(r_s, 0, len(edo_list) - 1):
                combo = [edo_list[i] for i in idx]
                t_mat: IntMat = np.vstack([k[1] for k in combo])
                t_mat = hnf(t_mat)

                t_tup = tuple(t_mat.flatten())
                if t_tup in checked:
                    continue

                checked.add(t_tup)

                label = " & ".join([k[0] for k in combo])

                if badness_type == "dirichlet":
                    b = temp_badness((t_mat, s), W_tenney_inv)
                else:
                    b = height(t_mat, W_weil_inv)

                by_rank[r_s].append((label, t_mat, b))
                total_count += 1
                count += 1
                if count >= 10:
                    break
            # print(count, math.comb(len(edo_list), r_s))

    res["results"] = ["families", "badness", "complexity"]
    for r_s in sorted(list(by_rank.keys())):
        res[f"= rank-{r_s} ="] = "<hr>"
        res_list = sorted(by_rank[r_s], key=lambda v: v[2])

        for k in res_list[0:10]:
            t_mat = np.atleast_2d(k[1])
            r, d = t_mat.shape
            compl = (d) ** (r - 1) * height(t_mat, W_tenney_inv) / np.sqrt(d)
            # compl = height(t_mat, W_tenney_inv) / np.sqrt(d)

            if s_non_prime:
                t_expanded = hnf(cokernel(basis @ kernel(t_mat)))
            else:
                t_expanded = t_mat

            families, families_weak = search_families(s_expanded, t_expanded)
            families_str = ""
            if len(families) > 0:
                families_str += ", ".join(sorted(list(families)))
            if len(families_weak) > 0 and len(families) == 0:
                if families_str != "":
                    families_str += ", "
                families_str += "(" + ", ".join(sorted(list(families_weak))) + ")"

            label = str(k[0])
            label_type = "edos"

            if r_s > dim - r_s:
                t_commas = kernel(t_mat)
                t_commas = LLL(t_commas, W=W_wilson)
                comma_rat = []
                for c in t_commas.T:
                    c = make_positive(c, s)
                    comma_str = str(ratio(c, s))
                    if len(comma_str) >= 11:
                        comma_str = "[" + " ".join(map(str, c)) + "]"
                    comma_rat.append(comma_str)
                label = ", ".join(comma_rat)
                label_type = "commas"

            url_args = {}
            if label_type == "edos":
                url_args["edos"] = label
                url_args["submit_edo"] = "submit"
            elif label_type == "commas":
                url_args["commas"] = label
                url_args["submit_comma"] = "submit"
            url_args["subgroup"] = args["subgroup"]
            url = Markup(f'<a href="{flask.url_for("result", **url_args)}">{label}</a>')

            res[url] = [families_str, f"{k[2]:.3f}", f"{compl:10.1f}"]

    print(f"search took: {time.time() - t_start} to find {total_count} temps, r={rank} d={dim}")

    return res
