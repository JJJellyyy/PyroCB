"""Actually RUN the matching pipeline that make_pcb_matching_clean_notebook.py writes
into PCB_MATCHING_CLEAN.ipynb.

make_pcb_matching_clean_notebook.py is a notebook GENERATOR (see its own docstring:
"Run: python make_pcb_matching_clean_notebook.py. Then open and run
PCB_MATCHING_CLEAN.ipynb in Jupyter") - its heavy computation lives as TEXT inside
`code(r"""...""")` calls, appended to `cells`, never executed by that script itself.
To regenerate pcb_matching_clean.gpkg/.tsv (e.g. after adding a new match version)
without opening Jupyter, this pulls every code cell's source out in the exact order
the notebook would run them and execs them in one shared namespace - so the notebook
source stays the single place the algorithm is written (no second copy to drift), and
running it here reproduces exactly what "run every cell top to bottom" would do.

Safe to re-run: relies on code() being called with a raw triple-double-quoted string
(the file's own convention throughout - see its module docstring), and on no cell
containing a literal \"\"\" of its own (function docstrings in that file use '''
specifically to avoid this collision - verified against the source, not assumed).
"""
import re
import time

SRC = "make_pcb_matching_clean_notebook.py"


def main():
    text = open(SRC, "r", encoding="utf-8").read()
    cells = re.findall(r'code\(r"""(.*?)"""\)', text, re.DOTALL)
    if not cells:
        raise RuntimeError(f"no code() cells found in {SRC} - has its format changed?")
    print(f"extracted {len(cells)} code cells from {SRC}", flush=True)

    ns: dict = {"__name__": "__pcb_matching_clean_run__"}
    t0 = time.time()
    for i, cell in enumerate(cells, 1):
        t1 = time.time()
        try:
            exec(compile(cell, f"<cell {i}>", "exec"), ns)
        except Exception:
            print(f"\n--- FAILED in cell {i}/{len(cells)} ---")
            raise
        print(f"  cell {i}/{len(cells)} ran in {time.time()-t1:5.1f}s", flush=True)
    print(f"\nall {len(cells)} cells ran in {time.time()-t0:.0f}s total")


if __name__ == "__main__":
    main()
