from pages import wiki_pages
import re

re_ratio = re.compile(r"^(\d+)/(\d+)$")

ratios: set[tuple[int, int]] = set()
for title in wiki_pages:
    m = re_ratio.match(title)
    if m is not None:
        p, q = m.groups()
        # print(f"{repr(title)} -> {p}/{q}")
        ratios.add((int(p), int(q)))

print("wiki_ratios = {")
ratios = sorted(list(ratios))
for k in ratios:
    print(f"\t{repr(k)},")
print("}")
