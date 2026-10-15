"""c6_sweep_tally_v1.py -- tally the C6 faithful-LP sweep (the 20 c6_lp_*.log
files written by cluster_cores_v4.py --lp-only i 20) and assert the counts quoted
in the paper: 384 cores; four-cores: 80 with LP value 1/3; five-cores: 192 with
LP value 1/4 and 112 with LP value 1/3.  Writes c6_above_cores.txt (the 112).
Usage: python c6_sweep_tally_v1.py <directory with c6_lp_*.log>"""
import sys, re, glob, os
from fractions import Fraction
from collections import Counter
rows = []
for f in glob.glob(os.path.join(sys.argv[1], 'c6_lp_*.log')):
    for line in open(f):
        m = re.match(r"idx (\d+) size (\d+) core (.*?): LP = (\S+) \(([\d.]+)\), exact-upper = (\S+), strategies = (\d+)", line)
        if m: rows.append((int(m[1]), int(m[2]), m[3], Fraction(m[4]), float(m[5])))
assert len(rows) == 384 and len(set(r[0] for r in rows)) == 384, len(rows)
assert max(abs(r[4] - float(r[3])) for r in rows) < 1e-6
c = Counter((r[1], str(r[3])) for r in rows)
print(dict(c))
assert c == {(4, '1/3'): 80, (5, '1/4'): 192, (5, '1/3'): 112}, c
above = sorted(r for r in rows if r[3] > Fraction(1, r[1] - 1))
with open('c6_above_cores.txt', 'w') as f:
    for r in above: f.write(f"{r[1]} {r[2].replace(' ', ',')}\n")
print("ASSERTED: 384 cores; 80 four-cores at 1/3; 192 five-cores at 1/4; 112 five-cores at 1/3 -> c6_above_cores.txt")
