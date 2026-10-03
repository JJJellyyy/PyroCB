"""Split fire_pcb_map.html into editable pieces and join them back.

    python map_tool.py split   -> map_src/head.html, map_src/data.json, map_src/tail.html
    python map_tool.py join    -> fire_pcb_map.html  (and prints the size)

The 47 MB map is one HTML file whose payload (`const DATA={...};`) is a single enormous line,
which makes it unreadable to editors and diff tools. head.html (CSS + markup) and tail.html
(all the JS) are the hand-written parts and are tracked in git; data.json is generated and is not.
"""
import json, sys

HTML, HEAD, DATA, TAIL = "fire_pcb_map.html", "map_src/head.html", "map_src/data.json", "map_src/tail.html"
MARK = "const DATA="

def split():
    t = open(HTML, encoding="utf-8", newline="").read()
    i = t.index(MARK)
    _, end = json.JSONDecoder().raw_decode(t[i + len(MARK):])
    j = i + len(MARK) + end
    open(HEAD, "w", encoding="utf-8", newline="\n").write(t[:i])
    open(DATA, "w", encoding="utf-8", newline="\n").write(t[i + len(MARK):j])
    open(TAIL, "w", encoding="utf-8", newline="\n").write(t[j:])
    print("head", i, "data", j - i - len(MARK), "tail", len(t) - j)

def join():
    out = open(HEAD, encoding="utf-8", newline="").read() + MARK + open(DATA, encoding="utf-8", newline="").read() + open(TAIL, encoding="utf-8", newline="").read()
    open(HTML, "w", encoding="utf-8", newline="\n").write(out)
    print("wrote", HTML, len(out))

if __name__ == "__main__":
    {"split": split, "join": join}[sys.argv[1]]()
