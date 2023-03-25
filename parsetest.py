import re
from lib_temper import *

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
ratio_pattern = r"(\d+)[/:](\d+)?"


subgroup = p_limit(11)

def parse_edos(s, subgroup):
    s = [int(i) for i in re.split(r"[\\.,; &]+", s.strip())]

    edos = []
    for e in s:
        edos.append(patent_map(e, subgroup))

    return edos

def parse_edos2(s, subgroup):
    # split string by separators (,.;& ) but ignore when in square brackets
    # "12, 17[+5, 8], 22"
    # => ["12", "17[+5, 8]", "22"]
    s = re.split(r"[\\.,; &]+(?!(?:[^,\[\]]+,)*[^,\[\]]+])", s.lower().strip())

    edos = []

    log_s = log_subgroup(subgroup)
    log_s = log_s / log_s[0]  # fix equave

    for e in s:
        print("parse: ", repr(e))
        res = re.split(r"(\d+)", e)[1:]
        if res == []:
            continue
        edo_num = int(res[0])
        p_map = patent_map(edo_num, subgroup)
        print(edo_num)
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
                    print(repr(ratio))

                    index = None
                    for i, pr in enumerate(subgroup):
                        if ratio == pr:
                            index = i
                            break
                    if index is None:
                        raise AssertionError("Adjustment not in subgroup")

                    for c in adj_str:
                        if c == "+":
                            p_map[0,index] = p_map[0,index] + 1
                        elif c == "-":
                            p_map[0,index] = p_map[0,index] - 1

                edos.append(p_map)
            else:
                # here we try to parse formats like '17c' (aka wart notation)
                warts = re.finditer(r"(\D)\1*", res[1])
                for g in warts:
                    w = g.group()
                    w_count = len(w)
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

                    p_map[0,index] = p_map[0,index] + sign * count

                edos.append(p_map)

    return edos

print(parse_edos2("156[-5, +++11] 156", subgroup))
# print(parse_edos("12 25 72", subgroup))


s = "12   522, 72aaabbce  & 156[-5, +++11]"




