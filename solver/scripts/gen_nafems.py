#!/usr/bin/env python3
"""Generate NAFEMS benchmark decks (LE1, LE10, T4) as CalculiX INP.

Lengths for LE1/LE10 are millimetres, stresses come out in MPa.
T4 is SI (metres, °C, W/m·K).

Targets (NAFEMS TNSB rev. 3 / Cameron, Casey & Simpson):
  LE1  σyy at D (inner, +x)           =  92.7 MPa
  LE10 σyy at D (inner, +x, top)      =  -5.38 MPa
  T4   T at (0.6 m, 0.2 m)            =  18.3 °C
"""
from __future__ import annotations

import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "examples" / "nafems"


def ellipse(s: float, th: float, ai, bi, ao, bo):
    a = ai + s * (ao - ai)
    b = bi + s * (bo - bi)
    return a * math.cos(th), b * math.sin(th)


def write_nodes(fh, ids, coords):
    fh.write("*NODE\n")
    for i in ids:
        x, y, z = coords[i]
        fh.write(f"{i}, {x:.8g}, {y:.8g}, {z:.8g}\n")


def write_set(fh, kind, name, members, per=12):
    fh.write(f"*{kind}, {kind}={name}\n")
    line = []
    for i, m in enumerate(members):
        line.append(str(m))
        if len(line) == per or i == len(members) - 1:
            fh.write(", ".join(line) + "\n")
            line = []


def gen_le1(nr=8, nt=24):
    """Quarter elliptic membrane, CPS8. Outward pressure 10 MPa on the outer edge."""
    ai, bi, ao, bo = 2000.0, 1000.0, 3250.0, 2750.0
    ni, nj = 2 * nr, 2 * nt
    coords = {}
    nid = {}
    n = 0
    for i in range(ni + 1):
        for j in range(nj + 1):
            if (i % 2) and (j % 2):
                continue  # face centre, not a serendipity node
            n += 1
            s = i / ni
            th = (j / nj) * math.pi / 2.0
            x, y = ellipse(s, th, ai, bi, ao, bo)
            nid[(i, j)] = n
            coords[n] = (x, y, 0.0)
    elems = []
    outer = []
    eid = 0
    for ie in range(nr):
        for je in range(nt):
            i0, j0 = 2 * ie, 2 * je
            conn = [
                nid[(i0, j0)],
                nid[(i0 + 2, j0)],
                nid[(i0 + 2, j0 + 2)],
                nid[(i0, j0 + 2)],
                nid[(i0 + 1, j0)],
                nid[(i0 + 2, j0 + 1)],
                nid[(i0 + 1, j0 + 2)],
                nid[(i0, j0 + 1)],
            ]
            eid += 1
            elems.append((eid, conn))
            if ie == nr - 1:
                outer.append(eid)
    d = nid[(0, 0)]
    path = OUT / "le1_cps8.inp"
    with path.open("w") as fh:
        fh.write("*HEADING\n")
        fh.write("NAFEMS LE1 elliptic membrane (plane stress, quarter)\n")
        fh.write("Target: syy(D) = 92.7 MPa at the inner point on +x\n")
        fh.write(f"** mesh CPS8 {nr} radial x {nt} hoop, node D = {d}\n")
        write_nodes(fh, range(1, n + 1), coords)
        fh.write("*ELEMENT, TYPE=CPS8, ELSET=MEM\n")
        for eid, c in elems:
            fh.write(str(eid) + ", " + ", ".join(map(str, c)) + "\n")
        write_set(fh, "ELSET", "OUTER", outer)
        write_set(fh, "NSET", "D", [d])
        ysym = [nid[(i, 0)] for i in range(ni + 1) if (i, 0) in nid]
        xsym = [nid[(i, nj)] for i in range(ni + 1) if (i, nj) in nid]
        write_set(fh, "NSET", "YS", ysym)
        write_set(fh, "NSET", "XS", xsym)
        fh.write("*MATERIAL, NAME=STEEL\n*ELASTIC\n210000, 0.3\n")
        fh.write("*SOLID SECTION, ELSET=MEM, MATERIAL=STEEL\n100\n")
        fh.write("*BOUNDARY\nYS, 2, 2\nXS, 1, 1\n")
        fh.write("*STEP\n*STATIC\n")
        # Positive P on a CCW edge is inward. Outward pressure → negative.
        fh.write("*DLOAD\nOUTER, P2, -10\n")
        fh.write("*NODE FILE\nU\n*EL FILE\nS\n*END STEP\n")
    return path, d, n, eid


def gen_le10(nr=4, nt=6, nz=2, name="le10_c3d20.inp"):
    """Thick elliptic plate, C3D20. NAFEMS fine mesh is 6 x 4 x 2 (hoop x radial x thick)."""
    ai, bi, ao, bo, h = 2000.0, 1000.0, 3250.0, 2750.0, 600.0
    ni, nj, nk = 2 * nr, 2 * nt, 2 * nz
    coords = {}
    nid = {}
    n = 0
    for i in range(ni + 1):
        for j in range(nj + 1):
            for k in range(nk + 1):
                odd = (i % 2) + (j % 2) + (k % 2)
                if odd > 1:
                    continue
                n += 1
                s = i / ni
                th = (j / nj) * math.pi / 2.0
                x, y = ellipse(s, th, ai, bi, ao, bo)
                z = (k / nk) * h
                nid[(i, j, k)] = n
                coords[n] = (x, y, z)
    elems = []
    top = []
    eid = 0
    for ie in range(nr):
        for je in range(nt):
            for ke in range(nz):
                i0, j0, k0 = 2 * ie, 2 * je, 2 * ke

                def g(di, dj, dk, i0=i0, j0=j0, k0=k0):
                    return nid[(i0 + di, j0 + dj, k0 + dk)]

                conn = [
                    g(0, 0, 0), g(2, 0, 0), g(2, 2, 0), g(0, 2, 0),
                    g(0, 0, 2), g(2, 0, 2), g(2, 2, 2), g(0, 2, 2),
                    g(1, 0, 0), g(2, 1, 0), g(1, 2, 0), g(0, 1, 0),
                    g(1, 0, 2), g(2, 1, 2), g(1, 2, 2), g(0, 1, 2),
                    g(0, 0, 1), g(2, 0, 1), g(2, 2, 1), g(0, 2, 1),
                ]
                eid += 1
                elems.append((eid, conn))
                if ke == nz - 1:
                    top.append(eid)
    d = nid[(0, 0, nk)]
    ysym = [nid[(i, 0, k)] for i in range(ni + 1) for k in range(nk + 1) if (i, 0, k) in nid]
    xsym = [nid[(i, nj, k)] for i in range(ni + 1) for k in range(nk + 1) if (i, nj, k) in nid]
    outer = [nid[(ni, j, k)] for j in range(nj + 1) for k in range(nk + 1) if (ni, j, k) in nid]
    omid = [nid[(ni, j, nz)] for j in range(nj + 1) if (ni, j, nz) in nid]
    path = OUT / name
    with path.open("w") as fh:
        fh.write("*HEADING\n")
        fh.write("NAFEMS LE10 thick plate pressure (quarter, C3D20)\n")
        fh.write("Target: syy(D) = -5.38 MPa at the inner top point on +x\n")
        fh.write(f"** mesh {nt} hoop x {nr} radial x {nz} thick, node D = {d}\n")
        write_nodes(fh, range(1, n + 1), coords)
        fh.write("*ELEMENT, TYPE=C3D20, ELSET=PLATE\n")
        for eid, c in elems:
            fh.write(str(eid) + ", " + ", ".join(map(str, c)) + "\n")
        write_set(fh, "ELSET", "TOP", top)
        write_set(fh, "NSET", "D", [d])
        write_set(fh, "NSET", "YS", ysym)
        write_set(fh, "NSET", "XS", xsym)
        write_set(fh, "NSET", "OUTER", outer)
        write_set(fh, "NSET", "OMID", omid)
        fh.write("*MATERIAL, NAME=STEEL\n*ELASTIC\n210000, 0.3\n")
        fh.write("*SOLID SECTION, ELSET=PLATE, MATERIAL=STEEL\n")
        fh.write("*BOUNDARY\n")
        fh.write("YS, 2, 2\nXS, 1, 1\nOUTER, 1, 2\nOMID, 3, 3\n")
        fh.write("*STEP\n*STATIC\n")
        # Positive P on the top face (P2) points outward (+z). Compression → negative.
        fh.write("*DLOAD\nTOP, P2, -1\n")
        fh.write("*NODE FILE\nU\n*EL FILE\nS\n*END STEP\n")
    return path, d, n, eid


def gen_t4(nx=24, ny=40, th=0.02):
    """NAFEMS thermal plate (Cameron/Casey/Simpson), thin C3D8 slice.

    0.6 m × 1.0 m, k=52 W/mK, bottom T=100°C, left insulated,
    top and right convection h=750 W/m²K to 0°C.
    Target T(0.6, 0.2) = 18.3°C.
    """
    # y-nodes include 0.2 exactly
    assert abs(round(0.2 / (1.0 / ny)) * (1.0 / ny) - 0.2) < 1e-12
    coords = {}
    nid = {}
    n = 0
    for ix in range(nx + 1):
        for iy in range(ny + 1):
            for iz in range(2):
                n += 1
                x = 0.6 * ix / nx
                y = 1.0 * iy / ny
                z = th * iz
                nid[(ix, iy, iz)] = n
                coords[n] = (x, y, z)
    right, top = [], []
    bottom = []
    probe = []
    eid = 0
    elems = []
    for ix in range(nx):
        for iy in range(ny):
            c = []
            for iz in (0, 1):
                for iy2, ix2 in ((iy, ix), (iy, ix + 1), (iy + 1, ix + 1), (iy + 1, ix)):
                    pass
            conn = [
                nid[(ix, iy, 0)],
                nid[(ix + 1, iy, 0)],
                nid[(ix + 1, iy + 1, 0)],
                nid[(ix, iy + 1, 0)],
                nid[(ix, iy, 1)],
                nid[(ix + 1, iy, 1)],
                nid[(ix + 1, iy + 1, 1)],
                nid[(ix, iy + 1, 1)],
            ]
            eid += 1
            elems.append((eid, conn))
            if ix == nx - 1:
                right.append(eid)
            if iy == ny - 1:
                top.append(eid)
    for iz in (0, 1):
        for ix in range(nx + 1):
            bottom.append(nid[(ix, 0, iz)])
        # probe on the right edge at y=0.2
        iy_p = int(round(0.2 * ny))
        probe.append(nid[(nx, iy_p, iz)])
    path = OUT / "t4_c3d8.inp"
    with path.open("w") as fh:
        fh.write("*HEADING\n")
        fh.write("NAFEMS thermal plate (T4 / Cameron, Casey & Simpson)\n")
        fh.write("Target: T(0.6 m, 0.2 m) = 18.3 C\n")
        fh.write(f"** mesh C3D8 {nx} x {ny} x 1, thickness {th} m, probe nodes {probe}\n")
        write_nodes(fh, range(1, n + 1), coords)
        fh.write("*ELEMENT, TYPE=C3D8, ELSET=PLATE\n")
        for eid, c in elems:
            fh.write(str(eid) + ", " + ", ".join(map(str, c)) + "\n")
        write_set(fh, "ELSET", "RIGHT", right)
        write_set(fh, "ELSET", "TOP", top)
        write_set(fh, "NSET", "BOTTOM", bottom)
        write_set(fh, "NSET", "PROBE", probe)
        fh.write("*MATERIAL, NAME=MAT\n*ELASTIC\n1, 0.3\n*CONDUCTIVITY\n52\n")
        fh.write("*SOLID SECTION, ELSET=PLATE, MATERIAL=MAT\n")
        fh.write("*BOUNDARY\nBOTTOM, 11, 11, 100\n")
        fh.write("*STEP\n*HEAT TRANSFER, STEADY STATE\n")
        fh.write("*FILM\nRIGHT, F4, 0, 750\nTOP, F5, 0, 750\n")
        fh.write("*NODE FILE\nNT\n*END STEP\n")
    return path, probe, n, eid


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for label, rec in (
        ("LE1", gen_le1()),
        ("LE10", gen_le10()),
        ("LE10f", gen_le10(8, 12, 4, "le10_c3d20_fine.inp")),
        ("T4", gen_t4()),
    ):
        path, probe, nnode, nelem = rec
        print(f"{label}: {path.name}  nodes={nnode}  elems={nelem}  probe={probe}")


if __name__ == "__main__":
    main()
