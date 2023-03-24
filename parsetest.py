import re

s = "12   522, 72aaabbce  & 156[-5, +++13]"

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

# split string by separators (,.;& ) but ignore when in square brackets
# "12, 17[+5, 8], 22"
# => ["12", "17[+5, 8]", "22"]
s = re.split(r"[\\.,; &]+(?!(?:[^,\[\]]+,)*[^,\[\]]+])", s.lower())

for e in s:
    print("parse: ", repr(e))
    res = re.split(r"(\d+)", e)[1:]
    edo_num = int(res[0])
    print(edo_num)
    if res[1] == "":
        pass
    else:
        adjust = re.findall(r"\[.*\]", e)
        if len(adjust) > 0:
            adjust = adjust[0][1:-1]
            adjust = re.split(r"[\\.,; &]+", adjust)
            for l in adjust:
                p = re.split(r"([+-]+)", l)[1:]
                print(p)
        else:
            warts = re.finditer(r"(\D)\1*", res[1])
            for g in (warts):
                print(g.group())
