from typing import Any
import time
import numpy as np
import flask

from info import from_commas, from_edos, search_families
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
    search_up = False
    if "edos" in args:
        t_map, _ = from_edos(args["edos"], basis, s_expanded)
        search_up = True
    elif "commas" in args:
        t_map, _ = from_commas(args["commas"], basis, s_expanded)
        search_up = False

    res: dict[str, Any] = {}

    rank = t_map.shape[0]
    dim = t_map.shape[1]

    if dim > 10:
        res["error"] = "please use a smaller subgroup"
        return res

    log_s = np.log2(np.array(s).astype(np.float64))

    ## metrics

    def subgroup_norm(n_mat: FloatMat) -> FloatMat:
        return (basis.T @ np.linalg.inv(n_mat) @ basis).astype(np.float64)

    def norm_det(n_mat: FloatMat) -> FloatMat:
        W_det = np.linalg.det(n_mat)
        assert W_det > 0
        return n_mat / np.power(W_det, 1 / n_mat.shape[0])

    # TE norm for badness
    W_tenney = subgroup_norm(metric_tenney(s_expanded))
    W_tenney_inv = np.linalg.inv(W_tenney)

    # Weil norm for cangwu badness
    W_weil = subgroup_norm(metric_weil_k(s_expanded, 800.0))
    W_weil_inv = np.linalg.inv(W_weil)

    # Wilson norm for displaying commas
    W_wilson = subgroup_norm(metric_wilson(s_expanded))
    W_wilson = W_wilson.astype(np.float64)

    # Weil norm for displaying edos
    W_edo = subgroup_norm(metric_weil_k(s_expanded, 1200.0 * np.max(log_s)))
    W_edo = np.linalg.inv(W_edo)

    # normalize determinants
    W_weil_inv = norm_det(W_weil_inv)
    W_tenney_inv = norm_det(W_tenney_inv)

    # pick some arbitrary edo to calibrate the complexity calculation
    compl_factor = 41.0 / height(patent_map(41.0, s), W_tenney_inv)

    t_start = time.time()
    total_count = 0

    # list of tuples (label, mapping, badness) ordered by rank
    by_rank: dict[int, list[tuple[str, IntMat, float]]] = {}

    checked: set[tuple[int, ...]] = set()

    def badness(t_map: IntMat) -> float:
        if badness_type == "dirichlet":
            return temp_badness((t_map, s), W_tenney_inv)
        return height(t_map, W_weil_inv)

    if search_up:
        if rank + 1 >= dim:
            res["warning"] = "empty search"
            return res

        for r_s in range(rank + 1, dim):
            by_rank[r_s] = []

        c_map = kernel(t_map)

        B = c_map

        f_init = 8 * np.max(log_s)
        f = 2 * (1.2 ** (dim - 1))

        # print(f"f init = {f_init}")
        # print(f"f = {f}")

        comma_list = []
        for i in range(8):
            try:
                W_LLL = subgroup_norm(metric_weil_k(s_expanded, f_init))
            except np.linalg.LinAlgError:
                break

            B = LLL(B, W=W_LLL, delta=0.9)
            f_init *= f
            for k in B.T:
                k = make_positive(k, s)  # TODO: speed up!

                label = str(ratio(k, s))
                if len(label) >= 11:
                    label = "[" + " ".join(map(str, k)) + "]"
                t_map = cokernel(k[:, None])  # can we calculate badness without this?

                b = badness(t_map)

                r_edo = (label, t_map, b)

                t_tup = tuple(k.flatten())
                if t_tup not in checked:
                    by_rank[dim - 1].append(r_edo)
                    comma_list.append((label, k, b))

                checked.add(t_tup)
        comma_list = sorted(comma_list, key=lambda v: v[2])

        max_per_rank = max(10, 170 // max(1, dim - 2))

        for r_s in range(2, dim - rank):
            count = 0
            for idx in comboBySum(r_s, 0, len(comma_list) - 1):
                combo = [comma_list[i] for i in idx]
                t_mat: IntMat = np.vstack([k[1] for k in combo]).T
                t_mat = cokernel(t_mat)
                if np.any(np.all(t_mat == 0, axis=0)):
                    continue

                t_tup = tuple(t_mat.flatten())
                if t_tup in checked:
                    continue

                checked.add(t_tup)

                label = ", ".join([k[0] for k in combo])

                b = badness(t_mat)

                by_rank[dim - r_s].append((label, t_mat, b))
                total_count += 1
                count += 1
                if count >= max_per_rank:
                    break

    else:
        if rank == 1:
            res["warning"] = "empty search"
            return res

        for r_s in range(1, rank):
            by_rank[r_s] = []

        f_init = 16 * np.sqrt(rank) * np.max(log_s)
        f = 1.4 * (1.1 ** (dim - 1))

        # print(f"f init = {f_init}")
        # print(f"f = {f}")

        B = t_map
        # f_prev = None

        edo_list = []
        for i in range(8):
            try:
                W_LLL = np.linalg.inv(subgroup_norm(metric_weil_k(s_expanded, f_init)))
            except np.linalg.LinAlgError:
                break
            B = LLL(B.T, W=W_LLL, delta=0.9).T
            f_init *= f

            found = np.abs(B)

            if np.any(found[:, 0] > 2000):
                break

            assert np.allclose(found - np.round(found), 0)
            found = np.round(found).astype(np.int64)

            # if f_prev and f_prev == tuple(found[:, 0]):
            #     break
            # f_prev = tuple(found[:, 0])
            # print(f_prev)

            for k in found:
                if 2 <= k[0] <= 1000:
                    t_edo = k[None, :]
                    label = edo_map_notation(k, s)

                    b = badness(k[None, :])

                    r_edo = (label, k, b)

                    t_tup = tuple(t_edo.flatten())
                    if t_tup not in checked:
                        by_rank[1].append(r_edo)
                        edo_list.append(r_edo)
                    checked.add(t_tup)

        max_per_rank = max(10, 170 // max(1, rank - 2))
        edo_list = sorted(edo_list, key=lambda v: v[2])
        # print(edo_list)

        for r_s in range(2, rank):
            count = 0
            # for combo in itertools.combinations(edo_list, r_s):
            for idx in comboBySum(r_s, 0, len(edo_list) - 1):
                combo = [edo_list[i] for i in idx]
                t_mat: IntMat = np.vstack([k[1] for k in combo])
                t_mat = hnf(t_mat)
                if np.any(np.all(t_mat[:] == 0, axis=1)):
                    continue

                t_tup = tuple(t_mat.flatten())
                if t_tup in checked:
                    continue

                checked.add(t_tup)

                label = " & ".join([k[0] for k in combo])

                b = badness(t_mat)

                by_rank[r_s].append((label, t_mat, b))
                total_count += 1
                count += 1
                if count >= max_per_rank:
                    break
            # print(count, math.comb(len(edo_list), r_s))

    res["results"] = ["families", "badness", "complexity"]
    for r_s in sorted(list(by_rank.keys())):
        res[f"= rank-{r_s} ="] = "<hr>"
        res_list = sorted(by_rank[r_s], key=lambda v: v[2])

        for k in res_list[0:10]:
            t_mat = np.atleast_2d(k[1])

            if factor_order(t_mat) > 1:
                continue

            r, d = t_mat.shape
            # compl = (d) ** (r - 1) * height(t_mat, W_tenney_inv) / np.sqrt(d)
            compl = (2) ** (r - 1) * height(t_mat, W_tenney_inv) * compl_factor

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

            if search_up:
                label_type = "commas"
                if r_s <= dim - r_s:
                    new_basis = LLL(t_mat.T, W=W_edo).T
                    labels = []
                    for edo_k in new_basis:
                        if edo_k[0] < 0:
                            edo_k *= -1
                        labels.append(edo_map_notation(edo_k, s))
                    label = " & ".join(labels)
                    label_type = "edos"

            else:
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

    s_str = ".".join([str(x) for x in s])
    print(f"search took: {time.time() - t_start} to find {total_count} temps, s={s_str}")

    return res
