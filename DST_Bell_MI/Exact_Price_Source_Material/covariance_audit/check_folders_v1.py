#!/usr/bin/env python3
"""
check_folders_v1.py -- which folder has which file, and does each match the public repo / the sent script?

Usage (Windows cmd):
  python check_folders_v1.py "P:\\CHSHDerivation\\moire-phase-space-sampler-main\\moire-phase-space-sampler-main\\DST_Bell_MI\\Exact_Price_Source_Material\\covariance_audit" "P:\\ApparentNonDeterminismLocal\\6\\breedrun"

Prints one row per file of interest: SHA-256 in each folder (or 'missing'), and whether it matches the
reference (the GitHub repo at commit a44817c, or the script as sent in chat).  Reference hashes were
computed on Oct 21, 2026.  Nothing is modified.
"""
import os, sys, hashlib

REF = {
    # repo inputs (commit a44817c, Sep 20 2026)
    "flight_time.py":                            "1641c4d93f82bc5bafd6913fe8274842b62c9244deebec18f58effe92ef7a707",
    "dither_sweep_rho.py":                       "2196890e23065c76625c5564f196d71e1fd27adb51d9650f7cfd23d56c402f61",
    "greedy_results.json":                       "e2d5368242393746c74a29f45888cad6761e1d659b6acb5bf35410beb358c673",
    "rho_nsh32_random0_eps0.01_nang161.json":    "e18b1ba7ca680fe74f644f8e2b2a94f914623bd9667af0245ba3d7cdc37c4434",
    "breed_t0.py":                               "0c4450a2a8ec35d7f1996176dbedb6dc4248d0d899eed2e9e5e52afb0ad1f93e",
    # scripts sent in chat, Oct 15-21
    "census_structure_v1.py":                    "e94eb4b4a59a22aca073319f2ea9d10d0389f5a4f128f680e8c76e165121794a",
    "census_structure_v2.py":                    "5c363e41ef215054d3fe9ce055142ddb87dbeb75ffed28ffa6cb2b3e1526d5b4",
    "terrain_factor_v1.py":                      "e1fb4ca877602e8e8604aab372093a3648144db81559225ec24fe9cdf635cf6f",
    "flight_time_perturb_v1.py":                 "47711d2b4e55b9c7234eb351bfef5ca35b4e813703b3a3eac236f6e33f9660be",
    "dither_sweep_perturb_v1.py":                "eaaad221bd5b09577802c42151322c63a1bdf69d0b027548da06706c86b89070",
    "product_census_v2.py":                      "d709a017ddc3bd6299f821862a1cedb99233f4bd5e23fed52a369468231a884e",
    "product_census_v3.py":                      "6ed535bcaa9d9947926dc10bd40813c20eaf1af354394ca40fff8ce35bc8bb02",
    "product_census_v4.py":                      "c6e7f2c3ee0f24fc97ef5fdfc55325771873f432eafe64141e8c87d1876393c0",
    "breed_product_v1.py":                       "b42b9a5886ea19cb1a7c754196c8b047bcef710bd76ed8335063c46e7f0ea65f",
    "pipeline_v1.py":                            "8955aa86b3d20e632f6f8aaa72d9799338497edcb46ba47e6adeef89656f1b1a",
}
# data files produced on Aaron's machine (no reference; just show which folder has them)
DATA = ["banks_nsh32_random0_nang161_nw2048.npz", "banks_nsh32_random1_nang161_nw2048.npz",
        "banks_nsh32_strat0_nang161_nw2048.npz", "banks_nsh32_strat1_nang161_nw2048.npz",
        "rho_nsh32_random1_eps0.01_nang161.json", "rho_nsh32_strat0_eps0.01_nang161.json",
        "rho_nsh32_strat1_eps0.01_nang161.json",
        "rho_perturb_2_2_1_shift_45_nsh32_random0_eps0.01_nang161.json",
        "rho_perturb_2_2_1_weight_0_nsh32_random0_eps0.01_nang161.json",
        "pipeline_FINAL.txt", "bred_product_champion.json"]

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()

def main(folders):
    names = [os.path.basename(os.path.normpath(f)) for f in folders]
    print("file".ljust(52) + "".join(n[:22].ljust(24) for n in names) + "verdict")
    for fn, ref in REF.items():
        row = fn.ljust(52); verdict = []
        for f, n in zip(folders, names):
            p = os.path.join(f, fn)
            if os.path.exists(p):
                h = sha(p); row += (h[:8] + ("  =ref" if h == ref else "  DIFF")).ljust(24)
                if h != ref: verdict.append(f"{n}: differs from reference")
            else:
                row += "missing".ljust(24)
        print(row + ("; ".join(verdict) if verdict else "ok"))
    print()
    print("data files (produced locally; no reference hash):")
    for fn in DATA:
        row = fn.ljust(52)
        for f in folders:
            p = os.path.join(f, fn)
            row += (("%.1f MB" % (os.path.getsize(p) / 1e6)) if os.path.exists(p) else "missing").ljust(24)
        print(row)
    # the decisive check: same greedy_results.json and flight_time.py in every folder that ran anything
    print()
    for f, n in zip(folders, names):
        for fn in ("greedy_results.json", "flight_time.py", "dither_sweep_rho.py"):
            p = os.path.join(f, fn)
            if os.path.exists(p) and sha(p) != REF[fn]:
                print(f"!! {n}\\{fn} is NOT the repo file -- any run in that folder used a different {fn}")

if __name__ == "__main__":
    main(sys.argv[1:])
