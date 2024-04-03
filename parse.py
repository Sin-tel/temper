# parsing and formatting utils
import re
import numpy as np
from typing import Optional
from markupsafe import Markup
from lib_temper import *
from wiki.ratios import wiki_ratios
from wiki.pages import wiki_pages


def cents(x, prec=3):
    return "{1:.{0}f}".format(prec, 1200 * x)


def ratio_with_link(frac: Fraction) -> str:
    frac_str = str(frac)
    p = (frac.numerator, frac.denominator)
    if p in wiki_ratios:
        return Markup(f'<a href="https://en.xen.wiki/w/{p[0]}/{p[1]}">{frac_str}</a>')
    return str(frac)


def page_with_link(titles: list[str], display: str) -> str:
    for title in titles:
        if title in wiki_pages:
            return Markup(f'<a href="https://en.xen.wiki/w/{title}">{display}</a>')
    return display


def parse_subgroup(s):
    s = [Fraction(i) for i in re.split(r"[\\.,; ]+", s)]

    # make sure they are all > 1
    for i, f in enumerate(s):
        p, q = f.as_integer_ratio()
        if np.log(p) - np.log(q) < 0.0:
            s[i] = Fraction(q, p)

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
            # if the input was simply an integer, then find the best mapping
            # TODO: this probably doesn't work correctly on subgroups
            best_b = 100000.0
            best_m = None
            search_range = (edo_num - 0.5, edo_num + 0.5)
            for m1 in Pmaps(search_range, subgroup):
                badness = temp_badness((m1, subgroup))
                if badness < best_b:
                    best_b = badness
                    best_m = np.copy(m1)
            if best_m is None:
                print(f"Somehow we did not find any patent maps for {edo_num} in {subgroup}")
                # fallback
                best_m = p_map
            edos.append(best_m)
        else:
            adjust = re.findall(r"\[.*\]", e)
            if len(adjust) > 0:
                # here we try to parse formats like '17[+5]'
                adjust = adjust[0][1:-1]
                adjust = re.split(r"[\\.,; &]+", adjust)
                for l in adjust:
                    p = re.split(r"([+-]+)", l)[1:]
                    adj_str = p[0]
                    rat = re.split(r"[/:]", p[1])
                    if len(rat) == 1:
                        rat = Fraction(int(rat[0]))
                    else:
                        rat = Fraction(int(rat[0]), int(rat[1]))

                    index = None
                    for i, pr in enumerate(subgroup):
                        if rat == pr:
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


# given a map for some edo, return the string representation in the new format
# (to replace warts notation)
def edo_map_notation(this_map, subgroup):
    this_edo = this_map[0]

    # divide by equave to fix non-octave temps
    j = log_subgroup(subgroup) / np.log2(float(subgroup[0]))

    p_map = np.round(this_edo * j).astype(np.int64)

    diff = this_map - p_map
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
