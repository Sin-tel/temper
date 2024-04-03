from typing import TypeAlias
import json
import numpy as np

from parse import parse_intervals
from lib_temper import *
from names import names

TempDict: TypeAlias = dict[tuple[int, ...], str]

# can be static if we don't generate it dynamically
# comma_filename = "temper/res/names.json"
comma_filename = "names.json"


def write_names() -> None:
    # can make this higher if needed
    s = p_limit(19)

    name_list: dict[str, dict[str, str]] = {}

    for name, comma_str in names:
        comma = parse_intervals(comma_str, s)
        assert len(comma) == 1
        comma = comma[0]
        # print(name)
        # print(ratio(comma, s))

        idx, _ = np.nonzero(comma)
        s_res = []
        for i in idx:
            s_res.append(s[i])
        c_res = comma[idx]
        # print(c_res)
        # print(s_res)
        s_key = str(tuple(s_res))

        if s_key not in name_list:
            name_list[s_key] = {}

        # already in hnf
        t = cokernel(c_res)
        # print(t)
        l_index = str(tuple(t.flatten()))
        name_list[s_key][l_index] = name
        # print("----")

    with open(comma_filename, "w", encoding="utf-8") as f:
        f.write(json.dumps(name_list, indent=4))


def parse_tuple(s_list: str) -> tuple[int, ...]:
    # parses tuple of integers, because json keys can't be tuples unlike in python
    #   "(1, 2, 3)" -> (1, 2, 3)
    return tuple(int(i) for i in s_list[1:-1].split(","))


def load_names() -> dict[tuple[int, ...], TempDict]:
    name_list: dict[tuple[int, ...], TempDict] = {}

    # convert json with string keys to tuple keys
    count = 0
    with open(comma_filename, "r", encoding="utf-8") as f:
        name_list_json = json.load(f)
        for k, v in name_list_json.items():
            key_subgroup = parse_tuple(k)
            d: TempDict = {}
            for t, name in v.items():
                key_t = parse_tuple(t)
                d[key_t] = name
                count += 1
            name_list[key_subgroup] = d

    print(f"loaded {count} temperament names!")
    return name_list


if __name__ == "__main__":
    write_names()
    load_names()
