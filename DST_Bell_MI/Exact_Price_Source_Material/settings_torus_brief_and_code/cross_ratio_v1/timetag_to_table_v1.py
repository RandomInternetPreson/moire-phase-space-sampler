#!/usr/bin/env python3
"""
timetag_to_table_v1.py — raw time-tags -> the coincidence table that
cross_ratio_v2.py (or v1) reads.  Companion to arXiv:2608.18886 (Sec. 10, Design A).
Python 3 + numpy.  Nothing here interprets the data; it only counts.

WHAT IT DOES
------------
Reads a time-tag stream (one record per detection: a time and a channel number),
removes a fixed delay per channel, and counts

  mode "table"  (four detectors per side): for each of the four basis
      combinations, the 2x2 coincidences between Alice's two detectors and
      Bob's two detectors in that basis, plus the singles on all eight
      detectors, the run duration and the coincidence window;

  mode "peaks"  (time-multiplexed basis choice, two detectors per user): for
      each pair of users, the 2x2 coincidences in the central peak (same basis)
      and in the two side peaks at +/- the long-path delay (cross basis), plus
      each user's two detector singles.

and writes the JSON that `python cross_ratio_v2.py out.json` analyses.

INPUT FORMATS (set in the config file)
--------------------------------------
  "text"    one detection per line, two numeric columns (time, channel) in any
            order, any delimiter; times in any unit ("time_unit_s").
  "binary"  fixed-length records described by a numpy dtype, e.g.
            [["time", "<u8"], ["chan", "<u2"]], after "header_bytes".
  "npz"     arrays "time" and "chan" in a numpy .npz file.
Files are read in chunks, so a file larger than memory is fine as long as it is
time-ordered.  Several files are taken as consecutive segments of one run
("sequential", the default) or as simultaneous streams on one clock
("parallel": loaded whole, merged and sorted).  Two taggers on separate clocks
must be synchronised first, as the lab's own sifting software already does.

USAGE
-----
  python timetag_to_table_v1.py config.json           # count, write the JSON
  python timetag_to_table_v1.py --scan config.json    # find per-channel delays
  python timetag_to_table_v1.py --example-config table
  python timetag_to_table_v1.py --example-config peaks
  python timetag_to_table_v1.py --selftest             # synthetic tags, both modes

Aaron Alai, November 2026.  ORCID 0000-0003-4704-7271.
"""
import json, os, sys, tempfile
import numpy as np

PS = 1e-12            # internal time unit: integer picoseconds

# ------------------------------------------------------------------ reading
def _to_ps(t, unit_s):
    """File time values -> int64 picoseconds."""
    t = np.asarray(t)
    if np.issubdtype(t.dtype, np.integer) and abs(unit_s / PS - round(unit_s / PS)) < 1e-9:
        return t.astype(np.int64) * int(round(unit_s / PS))
    return np.rint(t.astype(np.float64) * (unit_s / PS)).astype(np.int64)

def iter_chunks(spec, chunk_rows):
    """Yield (time_ps int64, chan int64) arrays from one input file, in order."""
    fmt = spec.get("format", "text")
    unit = float(spec.get("time_unit_s", 1e-12))
    path = spec["path"]
    if fmt == "text":
        tcol, ccol = int(spec.get("time_col", 0)), int(spec.get("chan_col", 1))
        delim = spec.get("delimiter", None)
        with open(path, "r") as fh:
            for _ in range(int(spec.get("skip_header", 0))):
                fh.readline()
            while True:
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    block = np.loadtxt(fh, delimiter=delim, usecols=(tcol, ccol),
                                       max_rows=chunk_rows, ndmin=2, dtype=np.float64)
                if block.shape[0] == 0:
                    break
                yield _to_ps(block[:, 0], unit), block[:, 1].astype(np.int64)
                if block.shape[0] < chunk_rows:
                    break
    elif fmt == "binary":
        dt = np.dtype([tuple(f) for f in spec["dtype"]])
        mm = np.memmap(path, dtype=dt, mode="r", offset=int(spec.get("header_bytes", 0)))
        for i in range(0, mm.shape[0], chunk_rows):
            blk = mm[i:i + chunk_rows]
            yield _to_ps(np.asarray(blk["time"]), unit), np.asarray(blk["chan"]).astype(np.int64)
    elif fmt == "npz":
        z = np.load(path)
        t, c = z["time"], z["chan"]
        for i in range(0, t.shape[0], chunk_rows):
            yield _to_ps(t[i:i + chunk_rows], unit), np.asarray(c[i:i + chunk_rows]).astype(np.int64)
    else:
        raise ValueError(f"unknown format {fmt!r}")

def stream(cfg):
    """Yield time-ordered (time_ps, chan) chunks for the whole run."""
    files = cfg["files"]
    rows = int(cfg.get("chunk_rows", 5_000_000))
    if cfg.get("merge", "sequential") == "parallel" and len(files) > 1:
        ts, cs = [], []
        for spec in files:
            off = int(round(float(spec.get("time_offset_s", 0.0)) / PS))
            choff = int(spec.get("chan_offset", 0))
            for t, c in iter_chunks(spec, rows):
                ts.append(t + off); cs.append(c + choff)
        t = np.concatenate(ts); c = np.concatenate(cs)
        order = np.argsort(t, kind="stable")
        yield t[order], c[order]
        return
    for spec in files:
        off = int(round(float(spec.get("time_offset_s", 0.0)) / PS))
        choff = int(spec.get("chan_offset", 0))
        for t, c in iter_chunks(spec, rows):
            yield t + off, c + choff

# ----------------------------------------------------------------- counting
def count_pairs(tA, tB, lo, hi, centre, half):
    """Number of (a, b) with lo <= a < hi and |b - a - centre| <= half.
    tA, tB sorted int64."""
    i0, i1 = np.searchsorted(tA, lo, "left"), np.searchsorted(tA, hi, "left")
    if i1 <= i0 or tB.size == 0:
        return 0
    a = tA[i0:i1] + centre
    return int(np.sum(np.searchsorted(tB, a + half, "right") - np.searchsorted(tB, a - half, "left")))

class Counter:
    """Streams chunks, keeps a tail so that pairs straddling a chunk boundary
    are counted exactly once, and accumulates the requested coincidences."""
    def __init__(self, channels, delays_ps, jobs, half_ps, margin_ps):
        self.channels = list(channels)              # channel numbers in use
        self.delays = {ch: int(delays_ps.get(ch, 0)) for ch in self.channels}
        self.jobs = jobs                            # list of (key, chanA, chanB, centre_ps)
        self.half = int(half_ps)
        self.margin = int(margin_ps)
        self.counts = {j[0]: 0 for j in jobs}
        self.singles = {ch: 0 for ch in self.channels}
        self.buf = {ch: np.empty(0, np.int64) for ch in self.channels}
        self.done_to = None                         # reference times < done_to are finished
        self.t_first = None; self.t_last = None

    def _process(self, upto):
        lo = self.done_to if self.done_to is not None else -(2 ** 62)
        for key, ca, cb, centre in self.jobs:
            self.counts[key] += count_pairs(self.buf[ca], self.buf[cb], lo, upto, centre, self.half)
        for ch in self.channels:
            b = self.buf[ch]
            self.singles[ch] += int(np.searchsorted(b, upto, "left") - np.searchsorted(b, lo, "left"))
        self.done_to = upto
        keep_from = upto - self.margin
        for ch in self.channels:
            b = self.buf[ch]
            self.buf[ch] = b[np.searchsorted(b, keep_from, "left"):]

    def feed(self, t, c):
        if t.size == 0:
            return
        if self.t_first is None:
            self.t_first = int(t[0])
        self.t_last = int(t[-1])
        for ch in self.channels:
            sel = t[c == ch] - self.delays[ch]
            if sel.size:
                if sel.size > 1 and np.any(np.diff(sel) < 0):
                    sel = np.sort(sel)
                self.buf[ch] = np.concatenate([self.buf[ch], sel])
        self._process(int(t[-1]) - self.margin)

    def finish(self):
        self._process(2 ** 62)
        dur = (self.t_last - self.t_first) * PS if self.t_first is not None else 0.0
        return dur

# -------------------------------------------------------------- mode: table
def run_table(cfg, quiet=False):
    ch = cfg["channels"]        # {"A_HV": [plus, minus], "A_DA": [...], "B_HV": [...], "B_DA": [...]}
    half = int(round(float(cfg["window_s"]) / PS / 2))
    delays = {int(k): int(round(float(v) / PS)) for k, v in cfg.get("delays_s", {}).items()}
    used = sorted({c for v in ch.values() for c in v})
    jobs = []
    for ab in ("HV", "DA"):
        for bb in ("HV", "DA"):
            for i in (0, 1):
                for j in (0, 1):
                    jobs.append(((f"{ab}_{bb}", i, j), ch["A_" + ab][i], ch["B_" + bb][j], 0))
    spread = (max(delays.values()) - min(delays.values())) if delays else 0
    C = Counter(used, delays, jobs, half, margin_ps=2 * half + spread + 1000)
    for t, c in stream(cfg):
        C.feed(t, c)
    dur = C.finish()
    out = {}
    for combo in ("HV_HV", "HV_DA", "DA_HV", "DA_DA"):
        out[combo] = [[C.counts[(combo, i, j)] for j in (0, 1)] for i in (0, 1)]
    out["singles"] = {k: [C.singles[v[0]], C.singles[v[1]]] for k, v in ch.items()}
    out["duration_s"] = dur
    out["window_s"] = 2 * half * PS
    if not quiet:
        print(f"table mode: {dur:.3f} s of data, window {2 * half} ps")
        for combo in ("HV_HV", "HV_DA", "DA_HV", "DA_DA"):
            print(f"  {combo}: {out[combo]}   sum {int(np.sum(out[combo]))}")
        print(f"  singles: {out['singles']}")
    return out

# -------------------------------------------------------------- mode: peaks
def run_peaks(cfg, quiet=False):
    users = cfg["users"]        # {"name": [plus_chan, minus_chan], ...}
    names = list(users)
    half = int(round(float(cfg["window_s"]) / PS / 2))
    Delta = int(round(float(cfg["long_path_delay_s"]) / PS))
    delays = {int(k): int(round(float(v) / PS)) for k, v in cfg.get("delays_s", {}).items()}
    pairs = cfg.get("pairs") or [[a, b] for i, a in enumerate(names) for b in names[i + 1:]]
    used = sorted({c for v in users.values() for c in v})
    jobs = []
    for a, b in pairs:
        for i in (0, 1):
            for j in (0, 1):
                # tau = t_b - t_a.  B on the long path arrives late: tau = +Delta.
                jobs.append(((a, b, "central", i, j), users[a][i], users[b][j], 0))
                jobs.append(((a, b, "side_ab", i, j), users[a][i], users[b][j], +Delta))
                jobs.append(((a, b, "side_ba", i, j), users[a][i], users[b][j], -Delta))
    spread = (max(delays.values()) - min(delays.values())) if delays else 0
    C = Counter(used, delays, jobs, half, margin_ps=2 * half + spread + abs(Delta) + 1000)
    for t, c in stream(cfg):
        C.feed(t, c)
    dur = C.finish()
    out = {"users": names, "pairs": [], "singles": {u: [C.singles[users[u][0]], C.singles[users[u][1]]] for u in names},
           "duration_s": dur, "window_s": 2 * half * PS}
    for a, b in pairs:
        rec = {"a": a, "b": b}
        for peak in ("central", "side_ab", "side_ba"):
            rec[peak] = [[C.counts[(a, b, peak, i, j)] for j in (0, 1)] for i in (0, 1)]
        out["pairs"].append(rec)
        if not quiet:
            print(f"  pair {a}-{b}: central {int(np.sum(rec['central']))}, "
                  f"side_ab (A short, B long) {int(np.sum(rec['side_ab']))}, "
                  f"side_ba (A long, B short) {int(np.sum(rec['side_ba']))}")
    if not quiet:
        print(f"peaks mode: {dur:.3f} s of data, window {2 * half} ps, long-path delay {Delta} ps")
    return out

# -------------------------------------------------------------------- scan
def scan(cfg):
    """Histogram t_B - t_A for every Alice/Bob channel pair on the first chunk
    and report where the coincidence peak sits, to set "delays_s"."""
    if cfg["mode"] == "table":
        A = sorted({c for k, v in cfg["channels"].items() if k.startswith("A_") for c in v})
        B = sorted({c for k, v in cfg["channels"].items() if k.startswith("B_") for c in v})
        pairs = [(a, b) for a in A for b in B]
    else:
        names = list(cfg["users"])
        pairs = [(ca, cb) for i, a in enumerate(names) for b in names[i + 1:]
                 for ca in cfg["users"][a] for cb in cfg["users"][b]]
    rng_ps = int(round(float(cfg.get("scan_range_s", 50e-9)) / PS))
    bin_ps = int(round(float(cfg.get("scan_bin_s", 100e-12)) / PS))
    t, c = next(stream(cfg))
    print(f"scan on the first {t.size} tags; range +/-{rng_ps} ps, bin {bin_ps} ps")
    print("  chanA chanB   peak at (ps)   peak counts   median bin")
    edges = np.arange(-rng_ps, rng_ps + bin_ps, bin_ps)
    for ca, cb in pairs:
        ta, tb = t[c == ca][:300000], t[c == cb]
        if ta.size == 0 or tb.size == 0:
            print(f"  {ca:5d} {cb:5d}   (no tags)"); continue
        lo = np.searchsorted(tb, ta - rng_ps, "left"); hi = np.searchsorted(tb, ta + rng_ps, "right")
        n = hi - lo
        idx = np.repeat(lo, n) + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
        d = tb[idx] - np.repeat(ta, n)
        h, _ = np.histogram(d, bins=edges)
        k = int(np.argmax(h))
        print(f"  {ca:5d} {cb:5d}   {int(edges[k] + bin_ps // 2):+10d}   {int(h[k]):11d}   {int(np.median(h)):10d}")
    print("Set delays_s so that every pair's peak sits at 0 (table) or at 0 and +/-long_path_delay (peaks):")
    print("a channel's delay is subtracted from its times, so delay_B - delay_A = peak position.")

# --------------------------------------------------------------- self-test
A_K = np.array([5.773047e-02, 9.834067e-03, 3.100629e-03, 1.331048e-03, 6.853601e-04, 3.973622e-04,
                2.503575e-04, 1.677139e-04, 1.177627e-04, 8.582274e-05, 6.445952e-05, 4.963526e-05])
def fringe(delta, share):
    k = np.arange(1, 13)
    return 1.0 + share * np.sum(A_K * np.cos(4 * k * delta))

def _synth_table_tags(rng, T_s, pair_rate, single_rate, share, delays_ps, jitter_ps, V=0.98, dark_rate=0.0):
    """Return (time_ps, chan, truth) for a four-detector-per-side link.
    Channels: A_HV 0,1  A_DA 2,3  B_HV 4,5  B_DA 6,7."""
    chan = {"A_HV": [0, 1], "A_DA": [2, 3], "B_HV": [4, 5], "B_DA": [6, 7]}
    eta = {"A_HV": [0.6, 0.4], "A_DA": [0.55, 0.45], "B_HV": [0.4, 0.6], "B_DA": [0.45, 0.55]}
    w = {"A": [0.55, 0.45], "B": [0.48, 0.52]}
    rel = {"HV_HV": 0.0, "DA_DA": 0.0, "HV_DA": np.pi / 4, "DA_HV": -np.pi / 4}
    sgn = np.array([[1, -1], [-1, 1]], float)
    ts, cs, truth = [], [], {}
    for combo, dl in rel.items():
        ab, bb = combo.split("_")
        E = -V * np.cos(2 * dl)
        wa = w["A"][0 if ab == "HV" else 1]; wb = w["B"][0 if bb == "HV" else 1]
        mean = pair_rate * T_s * wa * wb * 0.25 * np.outer(eta["A_" + ab], eta["B_" + bb]) * (1 + sgn * E) * fringe(dl, share)
        n = rng.poisson(mean); truth[combo] = n.tolist()
        for i in (0, 1):
            for j in (0, 1):
                t0 = rng.integers(0, int(T_s / PS), size=n[i, j])
                ca, cb = chan["A_" + ab][i], chan["B_" + bb][j]
                ts.append(t0 + delays_ps[ca] + np.rint(rng.normal(0, jitter_ps, n[i, j])).astype(np.int64)); cs.append(np.full(n[i, j], ca))
                ts.append(t0 + delays_ps[cb] + np.rint(rng.normal(0, jitter_ps, n[i, j])).astype(np.int64)); cs.append(np.full(n[i, j], cb))
    # unpaired detections (partner lost): proportional to splitter weight x efficiency;
    # dark counts: the same on every detector
    for key, pair in chan.items():
        wgt = w[key[0]][0 if key[2:] == "HV" else 1]
        for i, ch in enumerate(pair):
            n = rng.poisson((single_rate * 4.0 * wgt * eta[key][i] + dark_rate) * T_s)
            ts.append(rng.integers(0, int(T_s / PS), size=n)); cs.append(np.full(n, ch))
    t = np.concatenate(ts); c = np.concatenate(cs)
    order = np.argsort(t, kind="stable")
    return t[order], c[order].astype(np.int64), truth, chan

def _synth_peaks_tags(rng, T_s, pair_rate, single_rate, Delta_ps, delays_ps, jitter_ps, V=0.98):
    """Three users, two detectors each (channels 0..5), time-multiplexed basis."""
    users = {"A": [0, 1], "B": [2, 3], "C": [4, 5]}
    x = {"A": 0.8, "B": 0.9, "C": 0.7}
    eta = {"A": [0.52, 0.48], "B": [0.5, 0.5], "C": [0.47, 0.53]}
    sgn = np.array([[1, -1], [-1, 1]], float)
    names = list(users); ts, cs, truth = [], [], {}
    for ia, a in enumerate(names):
        for b in names[ia + 1:]:
            base = pair_rate * T_s * 0.25 * np.outer(eta[a], eta[b])
            cats = {"ss": (0, 0, 1.0, -V), "ll": (1, 1, x[a] * x[b], -V), "sl": (0, 1, x[b], 0.0), "ls": (1, 0, x[a], 0.0)}
            rec = {"central": np.zeros((2, 2), int), "side_ab": np.zeros((2, 2), int), "side_ba": np.zeros((2, 2), int)}
            for name, (la, lb, wgt, E) in cats.items():
                n = rng.poisson(base * wgt * (1 + sgn * E))
                key = "central" if la == lb else ("side_ab" if lb == 1 else "side_ba")
                rec[key] += n
                for i in (0, 1):
                    for j in (0, 1):
                        t0 = rng.integers(0, int(T_s / PS), size=n[i, j])
                        ca, cb = users[a][i], users[b][j]
                        ts.append(t0 + la * Delta_ps + delays_ps[ca] + np.rint(rng.normal(0, jitter_ps, n[i, j])).astype(np.int64)); cs.append(np.full(n[i, j], ca))
                        ts.append(t0 + lb * Delta_ps + delays_ps[cb] + np.rint(rng.normal(0, jitter_ps, n[i, j])).astype(np.int64)); cs.append(np.full(n[i, j], cb))
            truth[(a, b)] = {k: v.tolist() for k, v in rec.items()}
    for u, pair in users.items():
        for i, ch in enumerate(pair):
            n = rng.poisson(single_rate * 2.0 * eta[u][i] * T_s)
            ts.append(rng.integers(0, int(T_s / PS), size=n)); cs.append(np.full(n, ch))
    t = np.concatenate(ts); c = np.concatenate(cs)
    order = np.argsort(t, kind="stable")
    return t[order], c[order].astype(np.int64), truth, users

def _check(label, got, truth, capture, acc, tol_sigma=5.0):
    """Extracted count vs truth * capture + accidentals, in Poisson sigmas."""
    got = np.asarray(got, float); truth = np.asarray(truth, float)
    expect = truth * capture + acc
    z = (got - expect) / np.sqrt(np.maximum(expect, 1.0))
    worst = float(np.max(np.abs(z)))
    ok = worst < tol_sigma
    print(f"    {label:22s} extracted {int(got.sum()):9d}  expected {expect.sum():11.1f}  worst cell {worst:4.2f} sigma  {'ok' if ok else 'FAIL'}")
    return ok

def selftest():
    from math import erf, sqrt
    ok = True
    rng = np.random.default_rng(11)
    tmp = tempfile.mkdtemp(prefix="timetag_selftest_")
    print("=" * 78)
    print("SELF-TEST 1: table mode. 8 channels, 60 s, 5e4 pairs/s, unpaired singles ~1e4/s per")
    print("channel (proportional to efficiency) plus 2e3/s dark counts, per-channel delays up to")
    print("+/-5 ns, 150 ps jitter per detector, 1 ns window, planted share 0.10")
    delays = {ch: int(d) for ch, d in enumerate(rng.integers(-5000, 5000, size=8))}
    jitter = 150.0; T = 60.0
    t, c, truth, chan = _synth_table_tags(rng, T, 5.0e4, 1.0e4, 0.10, delays, jitter, dark_rate=2.0e3)
    print(f"  {t.size} tags generated")
    binpath = os.path.join(tmp, "tags.bin")
    rec = np.empty(t.size, dtype=[("time", "<u8"), ("chan", "<u2")]); rec["time"] = t + 10 ** 7; rec["chan"] = c; rec.tofile(binpath)
    base = {"mode": "table", "window_s": 1.0e-9, "channels": {k: v for k, v in chan.items()},
            "delays_s": {str(ch): d * PS for ch, d in delays.items()}}
    cfgA = dict(base, files=[{"path": binpath, "format": "binary", "dtype": [["time", "<u8"], ["chan", "<u2"]], "time_unit_s": 1e-12}], chunk_rows=10 ** 9)
    cfgB = dict(cfgA, chunk_rows=400_000)
    outA = run_table(cfgA, quiet=True); outB = run_table(cfgB, quiet=True)
    same = all(outA[k] == outB[k] for k in ("HV_HV", "HV_DA", "DA_HV", "DA_DA", "singles"))
    print(f"  one chunk vs {int(np.ceil(t.size / 400000))} chunks: tables identical -> {'ok' if same else 'FAIL'}"); ok &= same
    capture = erf(500.0 / (sqrt(2) * jitter * sqrt(2)))
    rate = {ch: outA["singles"][k][i] / T for k, v in chan.items() for i, ch in enumerate(v)}
    for combo in ("HV_HV", "HV_DA", "DA_HV", "DA_DA"):
        ab, bb = combo.split("_")
        acc = np.array([[rate[chan["A_" + ab][i]] * rate[chan["B_" + bb][j]] * 1e-9 * T for j in (0, 1)] for i in (0, 1)])
        ok &= _check(combo, outA[combo], truth[combo], capture, acc)
    R = {k: float(np.sum(outA[k])) for k in ("HV_HV", "HV_DA", "DA_HV", "DA_DA")}
    Rt = {k: float(np.sum(truth[k])) for k in R}
    X = R["HV_HV"] * R["DA_DA"] / (R["HV_DA"] * R["DA_HV"]); Xt = Rt["HV_HV"] * Rt["DA_DA"] / (Rt["HV_DA"] * Rt["DA_HV"])
    print(f"  raw cross-ratio from extracted table {X:.5f}; from the true pair list {Xt:.5f}")
    # text round trip on the first 5 s
    n5 = int(np.searchsorted(t, int(5.0 / PS)))
    txt = os.path.join(tmp, "tags.txt")
    np.savetxt(txt, np.column_stack([c[:n5], (t[:n5] + 10 ** 7) * 1e-12]), fmt=["%d", "%.12f"], header="channel time_s")
    cfgT = dict(base, files=[{"path": txt, "format": "text", "time_col": 1, "chan_col": 0, "time_unit_s": 1.0, "skip_header": 1}], chunk_rows=250_000)
    cfgU = dict(base, files=[{"path": binpath, "format": "binary", "dtype": [["time", "<u8"], ["chan", "<u2"]], "time_unit_s": 1e-12}], chunk_rows=10 ** 9)
    outT = run_table(cfgT, quiet=True)
    # the same 5 s from the binary file, via an npz slice
    npz = os.path.join(tmp, "tags5.npz"); np.savez(npz, time=t[:n5] + 10 ** 7, chan=c[:n5])
    outN = run_table(dict(base, files=[{"path": npz, "format": "npz", "time_unit_s": 1e-12}]), quiet=True)
    same = all(outT[k] == outN[k] for k in ("HV_HV", "HV_DA", "DA_HV", "DA_DA", "singles"))
    print(f"  text file (seconds, columns swapped, chunked) vs npz, first 5 s: identical -> {'ok' if same else 'FAIL'}"); ok &= same
    # delay scan recovers the planted delays
    print("  delay scan, first chunk (planted delay_B - delay_A in brackets):")
    t1, c1 = next(stream(dict(cfgU, chunk_rows=3_000_000)))
    worst = 0
    for ca, cb in ((0, 4), (2, 6), (1, 7), (3, 5)):
        ta, tb = t1[c1 == ca][:100000], t1[c1 == cb]
        lo = np.searchsorted(tb, ta - 20000, "left"); hi = np.searchsorted(tb, ta + 20000, "right"); n = hi - lo
        idx = np.repeat(lo, n) + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
        d = tb[idx] - np.repeat(ta, n)
        h, e = np.histogram(d, bins=np.arange(-20000, 20100, 100))
        pk = int(e[np.argmax(h)] + 50); want = delays[cb] - delays[ca]
        worst = max(worst, abs(pk - want))
        print(f"    channels {ca}-{cb}: peak at {pk:+6d} ps  [{want:+6d}]")
    good = worst <= 150
    print(f"    worst error {worst} ps (bin 100 ps) -> {'ok' if good else 'FAIL'}"); ok &= good

    print("=" * 78)
    print("SELF-TEST 2: peaks mode. 3 users x 2 detectors, long-path delay 3.7 ns, path ratios")
    print("0.8 / 0.9 / 0.7, 100 ps jitter, 1 ns window, 40 s, 3e4 pairs/s per user pair")
    delays = {ch: int(d) for ch, d in enumerate(rng.integers(-3000, 3000, size=6))}
    T = 40.0; jitter = 100.0
    t, c, truth, users = _synth_peaks_tags(rng, T, 3.0e4, 2.0e4, 3700, delays, jitter)
    print(f"  {t.size} tags generated")
    npz = os.path.join(tmp, "peaks.npz"); np.savez(npz, time=t + 10 ** 7, chan=c)
    cfgP = {"mode": "peaks", "window_s": 1.0e-9, "long_path_delay_s": 3.7e-9, "users": users,
            "delays_s": {str(ch): d * PS for ch, d in delays.items()},
            "files": [{"path": npz, "format": "npz", "time_unit_s": 1e-12}], "chunk_rows": 10 ** 9}
    outP = run_peaks(cfgP, quiet=True)
    outQ = run_peaks(dict(cfgP, chunk_rows=300_000), quiet=True)
    same = outP["pairs"] == outQ["pairs"] and outP["singles"] == outQ["singles"]
    print(f"  one chunk vs {int(np.ceil(t.size / 300000))} chunks: identical -> {'ok' if same else 'FAIL'}"); ok &= same
    capture = erf(500.0 / (sqrt(2) * jitter * sqrt(2)))
    rate = {ch: outP["singles"][u][i] / T for u, v in users.items() for i, ch in enumerate(v)}
    for recd in outP["pairs"]:
        a, b = recd["a"], recd["b"]
        acc = np.array([[rate[users[a][i]] * rate[users[b][j]] * 1e-9 * T for j in (0, 1)] for i in (0, 1)])
        for peak in ("central", "side_ab", "side_ba"):
            ok &= _check(f"{a}-{b} {peak}", recd[peak], truth[(a, b)][peak], capture, acc)
    # hand the two outputs to cross_ratio_v1.py if it sits next to this file
    here = os.path.dirname(os.path.abspath(__file__))
    cr = os.path.join(here, "cross_ratio_v2.py")
    if not os.path.exists(cr):
        cr = os.path.join(here, "cross_ratio_v1.py")
    for name, out in (("selftest_table_out.json", outA), ("selftest_peaks_out.json", outP)):
        with open(os.path.join(tmp, name), "w") as f:
            json.dump(out, f)
    if os.path.exists(cr):
        import importlib.util
        spec = importlib.util.spec_from_file_location("cross_ratio", cr)
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
        print("=" * 78); print(f"CHAIN: extracted table -> {os.path.basename(cr)} analyse_table (planted share 0.10, X = {m.X_of_share(0.10):.5f})")
        m.analyse_table(outA)
        print("-" * 78); print(f"CHAIN: extracted peaks -> {os.path.basename(cr)} analyse_peaks (planted share 0, path ratios 0.8/0.9/0.7)")
        m.analyse_peaks(outP)
    else:
        print("(cross_ratio_v2.py / v1 not found next to this file; chain test skipped)")
    print("=" * 78)
    print("SELF-TEST " + ("PASSED" if ok else "FAILED") + f"   (scratch files in {tmp})")
    return ok

EXAMPLE = {
"table": {"mode": "table", "window_s": 1.0e-9,
          "files": [{"path": "run01.txt", "format": "text", "time_col": 0, "chan_col": 1, "time_unit_s": 1e-12, "skip_header": 0}],
          "channels": {"A_HV": [1, 2], "A_DA": [3, 4], "B_HV": [5, 6], "B_DA": [7, 8]},
          "delays_s": {"1": 0.0, "2": 0.0, "3": 0.0, "4": 0.0, "5": 0.0, "6": 0.0, "7": 0.0, "8": 0.0},
          "output": "run01_table.json", "chunk_rows": 5000000},
"peaks": {"mode": "peaks", "window_s": 1.0e-9, "long_path_delay_s": 3.7e-9,
          "files": [{"path": "run01.bin", "format": "binary", "dtype": [["time", "<u8"], ["chan", "<u2"]], "header_bytes": 0, "time_unit_s": 1e-12}],
          "users": {"Alice": [1, 2], "Bob": [3, 4], "Chloe": [5, 6]},
          "delays_s": {"1": 0.0, "2": 0.0, "3": 0.0, "4": 0.0, "5": 0.0, "6": 0.0},
          "output": "run01_peaks.json", "chunk_rows": 5000000}}

if __name__ == "__main__":
    a = sys.argv[1:]
    if not a or a[0] in ("-h", "--help"):
        print(__doc__); sys.exit(0)
    if a[0] == "--selftest":
        sys.exit(0 if selftest() else 1)
    if a[0] == "--example-config":
        print(json.dumps(EXAMPLE[a[1]], indent=1)); sys.exit(0)
    do_scan = a[0] == "--scan"
    with open(a[1] if do_scan else a[0]) as f:
        cfg = json.load(f)
    if "channels" in cfg:
        cfg["channels"] = {k: [int(x) for x in v] for k, v in cfg["channels"].items()}
    if "users" in cfg:
        cfg["users"] = {k: [int(x) for x in v] for k, v in cfg["users"].items()}
    if do_scan:
        scan(cfg); sys.exit(0)
    out = run_table(cfg) if cfg["mode"] == "table" else run_peaks(cfg)
    dest = cfg.get("output", "timetag_table_out.json")
    with open(dest, "w") as f:
        json.dump(out, f, indent=1)
    print(f"wrote {dest}  ->  python cross_ratio_v2.py {dest}")
