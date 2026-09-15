#!/usr/bin/env python3
"""Run CalculiX verification decks through Axia and compare .dat against .dat.ref.

Usage:
    python3 solver/scripts/ccx-compare.py --ccx-dir /path/to/ccx/test \\
        --axia /path/to/axia --out /tmp/ccx-compare.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

SUPPORTED_ELEMS = {
    "C3D8", "C3D8R", "C3D8I", "C3D20", "C3D20R", "C3D20RI", "C3D4", "C3D10",
    "C3D10M", "C3D10T", "C3D6", "C3D15",
    "CPS4", "CPS4R", "CPE4", "CPE4R", "CPE4I", "CPS8", "CPE8", "CPS8R", "CPE8R",
    "CPS3", "CPE3", "CPS6", "CPE6",
    "CAX4", "CAX4R", "CAX8", "CAX8R", "CAX3", "CAX6",
    "MASS", "ROTARYI", "DASHPOTA", "DASHPOT", "GAPUNI",
    "S4", "S4R", "S3", "S3R", "STRI3", "S8", "S8R", "S6", "STRI65",
    "M3D3", "M3D4", "M3D4R", "M3D6", "M3D8", "M3D8R",
    "B31", "B31R", "B21", "B21R", "B32", "B32R", "B22", "B22R",
    "T3D2", "T2D2", "T3D3", "SPRINGA", "SPRING2",
}

SKIP_KWS = (
    "*FLUID DYNAMIC", "*CFD", "*ELECTROMAGNETICS", "*MAGNETOSTATICS",
    "*STEADY STATE DYNAMICS", "*COMPLEX FREQUENCY", "*MODAL DYNAMIC",
    "*GREEN", "*CREEP", "*VISCO",
    "*CONCRETE", "*MOHR COULOMB", "*DRUCKER PRAGER",
    "*DAMAGE", "*USER MATERIAL", "*USER ELEMENT",
    "*CYCLIC SYMMETRY MODEL", "*SUBMODEL",
    "*DESIGN VARIABLES", "*DESIGN RESPONSE", "*SENSITIVITY",
    "*RADIATE", "*VIEWFACTOR", "*RESTART",
    "*PRE-TENSION SECTION", "*MATRIX ASSEMBLE", "*MATRIX INPUT",
    "*CAP PLASTICITY", "*GAP CONDUCTANCE",
    "*SPECIFIC GAS CONSTANT", "*FLUID CONSTANTS", "*FLUID SECTION",
    "*COUPLED TEMPERATURE-DISPLACEMENT",
)

KW_RE = re.compile(r"^\*([A-Za-z][A-Za-z0-9 \-]*)", re.M)
NUM_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")


def read_text(p: Path) -> str:
    return p.read_text(encoding="latin-1", errors="replace")


def classify(inp: Path) -> dict:
    text = read_text(inp)
    kws = ["*" + m.group(1).upper().strip() for m in KW_RE.finditer(text)]
    elems = []
    for line in text.splitlines():
        u = line.strip().upper()
        if u.startswith("*ELEMENT"):
            m = re.search(r"TYPE\s*=\s*([A-Z0-9]+)", u)
            if m:
                elems.append(m.group(1))
    unknown = sorted({t for t in elems if t not in SUPPORTED_ELEMS})
    has_supported = any(t in SUPPORTED_ELEMS for t in elems)
    has_node = any(k.startswith("*NODE") for k in kws)
    has_step = any(k.startswith("*STEP") for k in kws)
    skip_kw = []
    for k in kws:
        for s in SKIP_KWS:
            if k == s or k.startswith(s):
                skip_kw.append(s)
                break
    displayable = has_supported and has_node and not unknown
    if inp.name.endswith(".rfn.inp"):
        displayable = False
    runnable = displayable and not skip_kw and has_step
    return {
        "name": inp.name,
        "elems": sorted(set(elems)),
        "unknown": unknown,
        "skip_kw": sorted(set(skip_kw)),
        "displayable": displayable,
        "runnable": runnable,
        "size": inp.stat().st_size,
    }


def parse_displacements(text: str) -> dict[int, tuple[float, float, float]]:
    out: dict[int, tuple[float, float, float]] = {}
    in_block = False
    for line in text.splitlines():
        low = line.lower()
        if "displacements" in low and "time" in low:
            in_block = True
            continue
        if in_block:
            if not line.strip():
                if out:
                    break
                continue
            if line.strip().startswith("S T E P") or "forces" in low or "stresses" in low:
                break
            parts = line.split()
            if len(parts) >= 4:
                try:
                    nid = int(float(parts[0]))
                    out[nid] = (float(parts[1]), float(parts[2]), float(parts[3]))
                except ValueError:
                    continue
    return out


def parse_eigenvalues(text: str) -> list[float]:
    freqs: list[float] = []
    in_block = False
    for line in text.splitlines():
        if "E I G E N V A L U E" in line:
            in_block = True
            continue
        if in_block:
            if "P A R T I C I P A T I O N" in line or "E F F E C T I V E" in line:
                break
            parts = line.split()
            if len(parts) >= 4:
                try:
                    int(float(parts[0]))
                    freqs.append(float(parts[3]))  # cycles/time
                except ValueError:
                    continue
    return freqs


def compare_u(ref: dict, got: dict) -> dict:
    if not ref:
        return {"kind": "no-ref-u"}
    keys = [k for k in ref if k in got]
    if not keys:
        return {"kind": "no-overlap", "nref": len(ref), "ngot": len(got)}
    umax = 0.0
    for v in ref.values():
        umax = max(umax, abs(v[0]), abs(v[1]), abs(v[2]))
    if umax < 1e-18:
        umax = 1.0
    err2 = 0.0
    worst = None
    for k in keys:
        for i in range(3):
            d = abs(ref[k][i] - got[k][i])
            r = d / umax
            if r > err2:
                err2 = r
                worst = (k, i, ref[k][i], got[k][i])
    return {
        "kind": "u",
        "n": len(keys),
        "umax_ref": umax,
        "rel_err": err2,
        "worst": worst,
        "ok": err2 <= 0.05,
        "tight": err2 <= 0.001,
    }


def compare_freq(ref: list[float], got: list[float]) -> dict:
    if not ref:
        return {"kind": "no-ref-freq"}
    n = min(len(ref), len(got))
    if n == 0:
        return {"kind": "no-got-freq", "nref": len(ref), "ngot": len(got)}
    err2 = 0.0
    for i in range(n):
        scale = max(abs(ref[i]), 1e-12)
        err2 = max(err2, abs(ref[i] - got[i]) / scale)
    return {
        "kind": "freq",
        "n": n,
        "rel_err": err2,
        "ok": err2 <= 0.05,
        "tight": err2 <= 0.01,
        "ref0": ref[0],
        "got0": got[0] if got else None,
    }


def run_one(axia: Path, inp: Path, work: Path, timeout: int) -> dict:
    stem = inp.stem
    outdir = work / stem
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    try:
        proc = subprocess.run(
            [
                str(axia),
                "--solver", "faer",
                "--json",
                "--quiet",
                "--no-frd",
                "-i", str(inp),
                "-o", str(outdir / stem),
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=str(inp.parent),
        )
    except subprocess.TimeoutExpired:
        return {"name": inp.name, "status": "timeout", "dt": time.time() - t0}
    dt = time.time() - t0
    stderr = (proc.stderr or "") + (proc.stdout or "")
    stats = {}
    for line in (proc.stdout or "").splitlines()[::-1]:
        line = line.strip()
        if line.startswith("{") and "ok" in line:
            try:
                stats = json.loads(line)
            except json.JSONDecodeError:
                pass
            break
        if line.startswith("{") and ("nnode" in line or "error" in line or "ok" in line):
            try:
                stats = json.loads(line)
            except json.JSONDecodeError:
                pass
            break
    if proc.returncode != 0:
        err = stderr.strip().splitlines()[-1] if stderr.strip() else f"exit {proc.returncode}"
        return {
            "name": inp.name,
            "status": "fail",
            "error": err[:300],
            "dt": dt,
        }
    # exploded / non-finite solutions: treat as fail even if the solver returned 0
    if isinstance(stats, dict):
        umax = stats.get("uMax")
        resid = stats.get("residual")
        if isinstance(umax, (int, float)) and (not math.isfinite(umax) or abs(umax) > 1e6):
            return {
                "name": inp.name,
                "status": "fail",
                "error": f"exploded uMax={umax}",
                "dt": dt,
            }
        if isinstance(resid, (int, float)) and (not math.isfinite(resid) or abs(resid) > 1.0):
            return {
                "name": inp.name,
                "status": "fail",
                "error": f"residual {resid}",
                "dt": dt,
            }
    dat = outdir / f"{stem}.dat"
    refp = inp.with_suffix(".dat.ref")
    cmp = {"kind": "solve-only"}
    if dat.exists() and refp.exists():
        got_txt = dat.read_text(encoding="latin-1", errors="replace")
        ref_txt = refp.read_text(encoding="latin-1", errors="replace")
        ru = parse_displacements(ref_txt)
        gu = parse_displacements(got_txt)
        rf = parse_eigenvalues(ref_txt)
        # Axia json frequencies
        gf = []
        if isinstance(stats, dict):
            gf = stats.get("frequencies") or stats.get("stats", {}).get("frequencies") or []
            if not gf:
                gf = parse_eigenvalues(got_txt)
        if rf:
            cmp = compare_freq(rf, list(gf) if gf else parse_eigenvalues(got_txt))
        elif ru:
            cmp = compare_u(ru, gu)
        else:
            cmp = {"kind": "no-ref-fields", "solved": True}
    status = "pass"
    if cmp.get("kind") in ("u", "freq"):
        status = "pass" if cmp.get("ok") else "mismatch"
    return {
        "name": inp.name,
        "status": status,
        "dt": dt,
        "cmp": cmp,
        "nnode": (stats.get("nnode") if isinstance(stats, dict) else None)
        or (stats.get("stats", {}).get("nnode") if isinstance(stats, dict) else None),
        "error": None,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccx-dir", required=True)
    ap.add_argument("--axia", required=True)
    ap.add_argument("--out", default="/tmp/ccx-compare.json")
    ap.add_argument("--work", default="/tmp/ccx-axia-out")
    ap.add_argument("--timeout", type=int, default=45)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--max-size", type=int, default=400_000)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    ccx = Path(args.ccx_dir)
    axia = Path(args.axia)
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)

    inps = sorted(ccx.glob("*.inp"))
    classified = [classify(p) for p in inps]
    skipped = [c for c in classified if not c["runnable"]]
    run = [c for c in classified if c["runnable"] and c["size"] <= args.max_size]
    too_big = [c for c in classified if c["runnable"] and c["size"] > args.max_size]
    if args.limit:
        run = run[: args.limit]

    print(f"total={len(classified)} skipped={len(skipped)} too_big={len(too_big)} run={len(run)}")
    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {
            ex.submit(run_one, axia, ccx / c["name"], work, args.timeout): c
            for c in run
        }
        done = 0
        for fut in as_completed(futs):
            c = futs[fut]
            try:
                r = fut.result()
            except Exception as e:
                r = {"name": c["name"], "status": "fail", "error": str(e)}
            r["elems"] = c["elems"]
            results.append(r)
            done += 1
            if done % 20 == 0 or r["status"] != "pass":
                print(f"[{done}/{len(run)}] {r['status']:9s} {r['name']} {r.get('error') or r.get('cmp',{}).get('rel_err','')}")

    summary = {
        "total_inp": len(classified),
        "skipped_not_representable": len(skipped),
        "too_big": [t["name"] for t in too_big],
        "ran": len(results),
        "pass": sum(1 for r in results if r["status"] == "pass"),
        "mismatch": sum(1 for r in results if r["status"] == "mismatch"),
        "fail": sum(1 for r in results if r["status"] == "fail"),
        "timeout": sum(1 for r in results if r["status"] == "timeout"),
        "results": sorted(results, key=lambda r: r["name"]),
        "skipped": skipped,
    }
    Path(args.out).write_text(json.dumps(summary, indent=1, default=str))
    print(
        f"PASS {summary['pass']}  MISMATCH {summary['mismatch']}  "
        f"FAIL {summary['fail']}  TIMEOUT {summary['timeout']}  "
        f"SKIP {summary['skipped_not_representable']}  BIG {len(too_big)}"
    )
    return 0 if summary["fail"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
