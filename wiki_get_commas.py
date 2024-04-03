from wiki_pages import wiki_pages
import re

re_ratio = re.compile(r"^(\d+)/(\d+)$")

ratios: list[tuple[int, int]] = []
for title in wiki_pages:
    m = re_ratio.match(title)
    if m is not None:
        p, q = m.groups()
        # print(f"{repr(title)} -> {p}/{q}")
        ratios.append((p, q))

print("wiki_ratios = ", repr(ratios))
