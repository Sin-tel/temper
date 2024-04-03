import urllib.request
import json
import time

headers = {
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Accept-Charset": "ISO-8859-1,utf-8;q=0.7,*;q=0.3",
    "Accept-Encoding": "none",
    "Accept-Language": "en-US,en;q=0.8",
    "Connection": "keep-alive",
}


def get(cont=None):
    limit = 500  # max is 500
    site = f"https://en.xen.wiki/api.php?action=query&list=allpages&aplimit={limit}&format=json"
    if cont is not None:
        site = f"https://en.xen.wiki/api.php?action=query&list=allpages&aplimit={limit}&format=json&apcontinue={cont}"

    request = urllib.request.Request(site, None, headers)
    page = urllib.request.urlopen(request).read()
    table = json.loads(page)
    return table


table = get()
titles: list[str] = []
while True:
    for p in table["query"]["allpages"]:
        title = p["title"]
        titles.append(title)

    if "continue" in table:
        cont = table["continue"]["apcontinue"]
        # print(cont)
        time.sleep(0.2)  # don't get banned
        table = get(cont)
    else:
        break

print(repr(titles))
