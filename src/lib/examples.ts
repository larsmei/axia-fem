export type Example = {
  id: string;
  name: string;
  blurb: string;
  inp: string;
};

function hexLattice(
  nx: number,
  ny: number,
  nz: number,
  lx: number,
  ly: number,
  lz: number,
): { nodes: string; elems: string; nnode: number; id: (i: number, j: number, k: number) => number } {
  const id = (i: number, j: number, k: number) =>
    1 + i + j * (nx + 1) + k * (nx + 1) * (ny + 1);
  let nodes = "";
  for (let k = 0; k <= nz; k++) {
    for (let j = 0; j <= ny; j++) {
      for (let i = 0; i <= nx; i++) {
        const x = (lx * i) / nx;
        const y = (ly * j) / ny;
        const z = (lz * k) / nz;
        nodes += `${id(i, j, k)}, ${x.toFixed(6)}, ${y.toFixed(6)}, ${z.toFixed(6)}\n`;
      }
    }
  }
  let elems = "";
  let e = 1;
  for (let k = 0; k < nz; k++) {
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const n0 = id(i, j, k);
        const n1 = id(i + 1, j, k);
        const n2 = id(i + 1, j + 1, k);
        const n3 = id(i, j + 1, k);
        const n4 = id(i, j, k + 1);
        const n5 = id(i + 1, j, k + 1);
        const n6 = id(i + 1, j + 1, k + 1);
        const n7 = id(i, j + 1, k + 1);
        elems += `${e}, ${n0}, ${n1}, ${n2}, ${n3}, ${n4}, ${n5}, ${n6}, ${n7}\n`;
        e += 1;
      }
    }
  }
  return { nodes, elems, nnode: (nx + 1) * (ny + 1) * (nz + 1), id };
}

function quadLattice(nx: number, ny: number, lx: number, ly: number) {
  const id = (i: number, j: number) => 1 + i + j * (nx + 1);
  let nodes = "";
  for (let j = 0; j <= ny; j++) {
    for (let i = 0; i <= nx; i++) {
      const x = (lx * i) / nx;
      const y = (ly * j) / ny;
      nodes += `${id(i, j)}, ${x.toFixed(6)}, ${y.toFixed(6)}, 0\n`;
    }
  }
  let elems = "";
  let e = 1;
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const n0 = id(i, j);
      const n1 = id(i + 1, j);
      const n2 = id(i + 1, j + 1);
      const n3 = id(i, j + 1);
      elems += `${e}, ${n0}, ${n1}, ${n2}, ${n3}\n`;
      e += 1;
    }
  }
  return { nodes, elems, id };
}

const tensionPatch = `*HEADING
Uniaxialer Zugstab — C3D8 Patch-Test
** Ein Hexaeder, E=210000, nu=0.3. Erwartung: ux(x=10)=0.01, Sxx=210.
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
5, 0, 0, 10
6, 10, 0, 10
7, 10, 10, 10
8, 0, 10, 10
*ELEMENT, TYPE=C3D8, ELSET=SOLID
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=SOLID, MATERIAL=STEEL
*NSET, NSET=FIXED
1, 4, 5, 8
*BOUNDARY
FIXED, 1, 1
1, 2, 3
4, 3, 3
*STEP
*STATIC
*CLOAD
2, 1, 5250
3, 1, 5250
6, 1, 5250
7, 1, 5250
*NODE FILE
U, RF
*EL FILE
S
*END STEP
`;

function tensionFine(): string {
  const nx = 6;
  const ny = 2;
  const nz = 2;
  const { nodes, elems, id } = hexLattice(nx, ny, nz, 60, 10, 10);
  const a = 10 * 10;
  const force = 210 * a;
  const nLoaded = (ny + 1) * (nz + 1);
  const fNode = force / nLoaded;
  let fix = "*NSET, NSET=FIXED\n";
  for (let k = 0; k <= nz; k++) {
    for (let j = 0; j <= ny; j++) fix += `${id(0, j, k)},\n`;
  }
  let load = "*NSET, NSET=LOADED\n";
  for (let k = 0; k <= nz; k++) {
    for (let j = 0; j <= ny; j++) load += `${id(nx, j, k)},\n`;
  }
  let cload = "*CLOAD\n";
  for (let k = 0; k <= nz; k++) {
    for (let j = 0; j <= ny; j++) cload += `${id(nx, j, k)}, 1, ${fNode.toFixed(6)}\n`;
  }
  return `*HEADING
Zugstab 3D — ${nx}×${ny}×${nz} C3D8, σ=210
*NODE
${nodes}*ELEMENT, TYPE=C3D8, ELSET=SOLID
${elems}*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=SOLID, MATERIAL=STEEL
${fix}${load}*BOUNDARY
FIXED, 1, 1
${id(0, 0, 0)}, 2, 3
${id(0, ny, 0)}, 3, 3
*STEP
*STATIC
${cload}*NODE FILE
U, RF
*EL FILE
S
*END STEP
`;
}

function cantilever2d(): string {
  const nx = 24;
  const ny = 6;
  const { nodes, elems, id } = quadLattice(nx, ny, 100, 10);
  const p = 100 / (ny + 1);
  let fix = "*NSET, NSET=FIX, GENERATE\n";
  fix += `1, ${id(0, ny)}, ${nx + 1}\n`;
  let cload = "*CLOAD\n";
  for (let j = 0; j <= ny; j++) cload += `${id(nx, j)}, 2, ${(-p).toFixed(6)}\n`;
  return `*HEADING
Kragträger 2D — CPS4, L=100, h=10, P=100
** Euler-Bernoulli: δ = PL³/(3EI) ≈ 1.905
*NODE
${nodes}*ELEMENT, TYPE=CPS4, ELSET=PLATE
${elems}*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=PLATE, MATERIAL=STEEL
1.0
${fix}*BOUNDARY
FIX, 1, 2
*STEP
*STATIC
${cload}*NODE FILE
U
*EL FILE
S
*END STEP
`;
}

function cantilever3d(): string {
  const nx = 12;
  const ny = 3;
  const nz = 3;
  const { nodes, elems, id } = hexLattice(nx, ny, nz, 100, 10, 10);
  const nTip = (ny + 1) * (nz + 1);
  const p = 100 / nTip;
  let fix = "*NSET, NSET=FIX\n";
  for (let k = 0; k <= nz; k++) {
    for (let j = 0; j <= ny; j++) fix += `${id(0, j, k)},\n`;
  }
  let cload = "*CLOAD\n";
  for (let k = 0; k <= nz; k++) {
    for (let j = 0; j <= ny; j++) cload += `${id(nx, j, k)}, 3, ${(-p).toFixed(6)}\n`;
  }
  return `*HEADING
Kragträger 3D — C3D8, Last in −Z
*NODE
${nodes}*ELEMENT, TYPE=C3D8, ELSET=SOLID
${elems}*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=SOLID, MATERIAL=STEEL
${fix}*BOUNDARY
FIX, 1, 3
*STEP
*STATIC
${cload}*NODE FILE
U, RF
*EL FILE
S
*END STEP
`;
}

function tetBlock(): string {
  const nx = 3;
  const ny = 3;
  const nz = 3;
  const { nodes, id } = hexLattice(nx, ny, nz, 12, 12, 12);
  // 6 tets per hex along diagonal 0-6
  let elems = "";
  let e = 1;
  const tets = (n: number[]) => {
    const [n0, n1, n2, n3, n4, n5, n6, n7] = n;
    const split = [
      [n0, n1, n2, n6],
      [n0, n2, n3, n6],
      [n0, n3, n7, n6],
      [n0, n7, n4, n6],
      [n0, n4, n5, n6],
      [n0, n5, n1, n6],
    ];
    for (const t of split) {
      elems += `${e}, ${t[0]}, ${t[1]}, ${t[2]}, ${t[3]}\n`;
      e += 1;
    }
  };
  for (let k = 0; k < nz; k++) {
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        tets([
          id(i, j, k),
          id(i + 1, j, k),
          id(i + 1, j + 1, k),
          id(i, j + 1, k),
          id(i, j, k + 1),
          id(i + 1, j, k + 1),
          id(i + 1, j + 1, k + 1),
          id(i, j + 1, k + 1),
        ]);
      }
    }
  }
  const a = 12 * 12;
  const force = 210 * a;
  const nLoaded = (ny + 1) * (nz + 1);
  const fNode = force / nLoaded;
  let fix = "*NSET, NSET=FIXED\n";
  for (let k = 0; k <= nz; k++) for (let j = 0; j <= ny; j++) fix += `${id(0, j, k)},\n`;
  let cload = "*CLOAD\n";
  for (let k = 0; k <= nz; k++) {
    for (let j = 0; j <= ny; j++) cload += `${id(nx, j, k)}, 1, ${fNode.toFixed(6)}\n`;
  }
  return `*HEADING
Quader aus C3D4 — Zug in X, σ=210
*NODE
${nodes}*ELEMENT, TYPE=C3D4, ELSET=SOLID
${elems}*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=SOLID, MATERIAL=STEEL
${fix}*BOUNDARY
FIXED, 1, 1
${id(0, 0, 0)}, 2, 3
${id(0, ny, 0)}, 3, 3
*STEP
*STATIC
${cload}*NODE FILE
U, RF
*EL FILE
S
*END STEP
`;
}

function plateWithHole(): string {
  const r0 = 10;
  const w = 40;
  const nr = 8;
  const nt = 10;
  const id = (i: number, j: number) => 1 + i + j * (nr + 1);
  let nodes = "";
  for (let j = 0; j <= nt; j++) {
    const theta = (Math.PI / 2) * (j / nt);
    const c = Math.cos(theta);
    const s = Math.sin(theta);
    const ix = r0 * c;
    const iy = r0 * s;
    let ox: number;
    let oy: number;
    if (theta <= Math.PI / 4) {
      ox = w;
      oy = w * Math.tan(theta);
    } else {
      ox = w * Math.tan(Math.PI / 2 - theta);
      oy = w;
    }
    for (let i = 0; i <= nr; i++) {
      const t = i / nr;
      // cluster near the hole
      const g = t * t;
      const x = ix + (ox - ix) * g;
      const y = iy + (oy - iy) * g;
      nodes += `${id(i, j)}, ${x.toFixed(6)}, ${y.toFixed(6)}, 0\n`;
    }
  }
  let elems = "";
  let e = 1;
  for (let j = 0; j < nt; j++) {
    for (let i = 0; i < nr; i++) {
      const n0 = id(i, j);
      const n1 = id(i + 1, j);
      const n2 = id(i + 1, j + 1);
      const n3 = id(i, j + 1);
      elems += `${e}, ${n0}, ${n1}, ${n2}, ${n3}\n`;
      e += 1;
    }
  }
  let symX = "*NSET, NSET=SYMX\n";
  for (let i = 0; i <= nr; i++) symX += `${id(i, 0)},\n`;
  let symY = "*NSET, NSET=SYMY\n";
  for (let i = 0; i <= nr; i++) symY += `${id(i, nt)},\n`;
  let far = "*NSET, NSET=FAR\n";
  for (let j = 0; j <= nt; j++) far += `${id(nr, j)},\n`;
  // tension in x on the right edge (theta near 0, outer)
  // FAR includes top of square too. Apply ux or CLOAD only on x=w edge.
  let cload = "*CLOAD\n";
  const pressure = 10;
  // approximate nodal forces on x=w vertical edge (j from 0 until y < w-eps)
  const edge: number[] = [];
  for (let j = 0; j <= nt; j++) {
    const theta = (Math.PI / 2) * (j / nt);
    if (theta <= Math.PI / 4 + 1e-9) edge.push(id(nr, j));
  }
  const f = (pressure * w) / Math.max(1, edge.length);
  for (const n of edge) cload += `${n}, 1, ${f.toFixed(6)}\n`;
  return `*HEADING
Viertelplatte mit Loch — CPS4, σ∞=10, Kt≈3
*NODE
${nodes}*ELEMENT, TYPE=CPS4, ELSET=PLATE
${elems}*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=PLATE, MATERIAL=STEEL
1.0
${symX}${symY}${far}*BOUNDARY
SYMX, 2, 2
SYMY, 1, 1
*STEP
*STATIC
${cload}*NODE FILE
U
*EL FILE
S
*END STEP
`;
}

function cantileverB32(): string {
  const l = 200;
  const nel = 8;
  const nnode = 2 * nel + 1;
  let nodes = "";
  for (let i = 0; i < nnode; i++) {
    const x = (l * i) / (nnode - 1);
    nodes += `${i + 1}, ${x.toFixed(4)}, 0, 0\n`;
  }
  let elems = "";
  for (let e = 0; e < nel; e++) {
    const n1 = 2 * e + 1;
    const mid = n1 + 1;
    const n2 = n1 + 2;
    elems += `${e + 1}, ${n1}, ${n2}, ${mid}\n`;
  }
  const tip = nnode;
  return `*HEADING
Kragträger B32 — Timoshenko, wie Abaqus
** L=200, RECT 10×20 (n1=Z, Höhe Y), E=210000, P=100 in −Y
** δ = PL³/3EI + PL/kAG ≈ 0.1906
*NODE
${nodes}*ELEMENT, TYPE=B32, ELSET=BEAM
${elems}*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*BEAM SECTION, ELSET=BEAM, MATERIAL=STEEL, SECTION=RECT
10, 20
0, 0, 1
*BOUNDARY
1, 1, 6
*STEP
*STATIC
*CLOAD
${tip}, 2, -100
*NODE FILE
U, RF
*EL FILE
S
*END STEP
`;
}

function portalB32(): string {
  const segs = (id0: number, e0: number, x0: number, y0: number, x1: number, y1: number) => {
    const nodes: string[] = [];
    const elems: string[] = [];
    const n = 5;
    for (let i = 0; i < n; i++) {
      const t = i / (n - 1);
      nodes.push(
        `${id0 + i}, ${(x0 + (x1 - x0) * t).toFixed(4)}, ${(y0 + (y1 - y0) * t).toFixed(4)}, 0`,
      );
    }
    elems.push(`${e0}, ${id0}, ${id0 + 2}, ${id0 + 1}`);
    elems.push(`${e0 + 1}, ${id0 + 2}, ${id0 + 4}, ${id0 + 3}`);
    return { nodes, elems };
  };
  const c1 = segs(1, 1, 0, 0, 0, 200);
  const c2 = segs(6, 3, 400, 0, 400, 200);
  return `*HEADING
Portalrahmen B32 — Stützen + Riegel, Last +X am Kopf
*NODE
${c1.nodes.join("\n")}
${c2.nodes.join("\n")}
11, 100, 200, 0
12, 200, 200, 0
13, 300, 200, 0
*ELEMENT, TYPE=B32, ELSET=FRAME
${c1.elems.join("\n")}
${c2.elems.join("\n")}
5, 5, 12, 11
6, 12, 10, 13
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*BEAM SECTION, ELSET=FRAME, MATERIAL=STEEL, SECTION=RECT
8, 16
0, 0, 1
*NSET, NSET=FIX
1, 6
*BOUNDARY
FIX, 1, 6
*STEP
*STATIC
*CLOAD
10, 1, 50
*NODE FILE
U, RF
*EL FILE
S
*END STEP
`;
}

function cantileverCps8(): string {
  const nx = 8;
  const ny = 2;
  const lx = 100;
  const ly = 10;
  let nodes = "";
  let next = 1;
  const idOf: number[][] = Array.from({ length: 2 * ny + 1 }, () =>
    Array<number>(2 * nx + 1).fill(0),
  );
  for (let j = 0; j <= 2 * ny; j++) {
    for (let i = 0; i <= 2 * nx; i++) {
      if (i % 2 === 1 && j % 2 === 1) continue;
      const x = (lx * i) / (2 * nx);
      const y = (ly * j) / (2 * ny);
      idOf[j][i] = next;
      nodes += `${next}, ${x.toFixed(6)}, ${y.toFixed(6)}, 0\n`;
      next += 1;
    }
  }
  let elems = "";
  let e = 1;
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const i0 = 2 * i;
      const j0 = 2 * j;
      const n1 = idOf[j0][i0];
      const n2 = idOf[j0][i0 + 2];
      const n3 = idOf[j0 + 2][i0 + 2];
      const n4 = idOf[j0 + 2][i0];
      const n5 = idOf[j0][i0 + 1];
      const n6 = idOf[j0 + 1][i0 + 2];
      const n7 = idOf[j0 + 2][i0 + 1];
      const n8 = idOf[j0 + 1][i0];
      elems += `${e}, ${n1}, ${n2}, ${n3}, ${n4}, ${n5}, ${n6}, ${n7}, ${n8}\n`;
      e += 1;
    }
  }
  let fix = "*NSET, NSET=FIX\n";
  for (let j = 0; j <= 2 * ny; j++) if (idOf[j][0]) fix += `${idOf[j][0]},\n`;
  const tip: number[] = [];
  for (let j = 0; j <= 2 * ny; j++) if (idOf[j][2 * nx]) tip.push(idOf[j][2 * nx]);
  const p = 100 / tip.length;
  let cload = "*CLOAD\n";
  for (const id of tip) cload += `${id}, 2, ${(-p).toFixed(6)}\n`;
  return `*HEADING
Kragträger CPS8 — quadratische Vierecke
** L=100, h=10, t=1, P=100. Euler δ ≈ 1.905
*NODE
${nodes}*ELEMENT, TYPE=CPS8, ELSET=PLATE
${elems}*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=PLATE, MATERIAL=STEEL
1.0
${fix}*BOUNDARY
FIX, 1, 2
*STEP
*STATIC
${cload}*NODE FILE
U
*EL FILE
S
*END STEP
`;
}

const tensionC3d20 = `*HEADING
Uniaxialer Zugstab — C3D20 Patch-Test
** Ein quadratisches Hexaeder, 20 Knoten. Erwartung: ux=0.01, Sxx=210.
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
5, 0, 0, 10
6, 10, 0, 10
7, 10, 10, 10
8, 0, 10, 10
9, 5, 0, 0
10, 10, 5, 0
11, 5, 10, 0
12, 0, 5, 0
13, 5, 0, 10
14, 10, 5, 10
15, 5, 10, 10
16, 0, 5, 10
17, 0, 0, 5
18, 10, 0, 5
19, 10, 10, 5
20, 0, 10, 5
*ELEMENT, TYPE=C3D20, ELSET=SOLID
1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=SOLID, MATERIAL=STEEL
*NSET, NSET=FIXED
1, 4, 5, 8, 12, 16, 17, 20
*BOUNDARY
FIXED, 1, 1
1, 2, 3
4, 3, 3
*STEP
*STATIC
*DLOAD
1, P4, 210
*NODE FILE
U, RF
*EL FILE
S
*END STEP
`;

function cantileverS4(): string {
  const nx = 8;
  const ny = 2;
  const { nodes, elems, id } = quadLattice(nx, ny, 100, 10);
  let fix = "";
  for (let j = 0; j <= ny; j++) fix += `${id(0, j)},\n`;
  const p = (1 / (ny + 1)).toFixed(6);
  let cload = "";
  for (let j = 0; j <= ny; j++) cload += `${id(nx, j)}, 3, ${-Number(p)}\n`;
  return `*HEADING
Kragplatte S4R — Mindlin/MITC4, 6 DOF
** L=100, b=10, t=1, E=210000, ν=0, P=1. Euler δ = PL³/3EI ≈ 1.905
*NODE
${nodes}*ELEMENT, TYPE=S4R, ELSET=PLATE
${elems}*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.0
*SHELL SECTION, ELSET=PLATE, MATERIAL=STEEL
1.0
*NSET, NSET=FIX
${fix}*BOUNDARY
FIX, 1, 6
*STEP
*STATIC
*CLOAD
${cload}*NODE FILE
U
*EL FILE
S
*END STEP
`;
}

function ssPlateS4(): string {
  const n = 8;
  const a = 100;
  const { nodes, elems, id } = quadLattice(n, n, a, a);
  // simply supported: w=0 on all edges, plus in-plane pin
  let bc = "";
  for (let i = 0; i <= n; i++) {
    bc += `${id(i, 0)}, 3, 3\n`;
    bc += `${id(i, n)}, 3, 3\n`;
  }
  for (let j = 1; j < n; j++) {
    bc += `${id(0, j)}, 3, 3\n`;
    bc += `${id(n, j)}, 3, 3\n`;
  }
  bc += `${id(0, 0)}, 1, 2\n`;
  bc += `${id(n, 0)}, 2, 2\n`;
  return `*HEADING
Gelagerte Quadratplatte S4R — Gleichlast
** a=100, t=1, q=0.01, E=210000, ν=0.3
** Kirchhoff δ_max ≈ 0.00406 q a⁴ / D ≈ 0.211
*NODE
${nodes}*ELEMENT, TYPE=S4R, ELSET=PLATE
${elems}*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SHELL SECTION, ELSET=PLATE, MATERIAL=STEEL
1.0
*BOUNDARY
${bc}*STEP
*STATIC
*DLOAD
EALL, P, -0.01
*NODE FILE
U
*EL FILE
S
*END STEP
`;
}

export const EXAMPLES: Example[] = [
  {
    id: "b32",
    name: "Kragträger B32",
    blurb: "Timoshenko-Balken, Abaqus B32",
    inp: cantileverB32(),
  },
  {
    id: "s4r",
    name: "Kragplatte S4R",
    blurb: "Mindlin-Schale MITC4, 6 DOF",
    inp: cantileverS4(),
  },
  {
    id: "plate",
    name: "Quadratplatte S4R",
    blurb: "Allseitig gelagert, Gleichlast",
    inp: ssPlateS4(),
  },
  {
    id: "frame",
    name: "Rahmen B32",
    blurb: "Portalrahmen aus B32",
    inp: portalB32(),
  },
  {
    id: "patch",
    name: "Zugstab C3D8",
    blurb: "Patch-Test, ein Hexaeder",
    inp: tensionPatch.trimStart(),
  },
  {
    id: "c3d15",
    name: "Wedge C3D15",
    blurb: "quadratischer Pentaeder, confined εz",
    inp: `*HEADING
C3D15 confined uniaxial strain
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 0, 10, 0
4, 0, 0, 10
5, 10, 0, 10
6, 0, 10, 10
7, 5, 0, 0
8, 5, 5, 0
9, 0, 5, 0
10, 5, 0, 10
11, 5, 5, 10
12, 0, 5, 10
13, 0, 0, 5
14, 10, 0, 5
15, 0, 10, 5
*ELEMENT, TYPE=C3D15, ELSET=S
1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=S, MATERIAL=STEEL
*BOUNDARY
1, 1, 2
2, 1, 2
3, 1, 2
4, 1, 2
5, 1, 2
6, 1, 2
7, 1, 2
8, 1, 2
9, 1, 2
10, 1, 2
11, 1, 2
12, 1, 2
13, 1, 2
14, 1, 2
15, 1, 2
1, 3, 3, 0.0
2, 3, 3, 0.0
3, 3, 3, 0.0
7, 3, 3, 0.0
8, 3, 3, 0.0
9, 3, 3, 0.0
13, 3, 3, 0.005
14, 3, 3, 0.005
15, 3, 3, 0.005
4, 3, 3, 0.01
5, 3, 3, 0.01
6, 3, 3, 0.01
10, 3, 3, 0.01
11, 3, 3, 0.01
12, 3, 3, 0.01
*STEP
*STATIC
*NODE FILE
U
*EL FILE
S
*END STEP
`,
  },
  {
    id: "cax4",
    name: "Rohr CAX4",
    blurb: "Achsensymmetrie, Innendruck, Lamé",
    inp: `*HEADING
CAX4 thick cylinder
*NODE
1, 10, 0
2, 10, 2
3, 12.5, 0
4, 12.5, 2
5, 15, 0
6, 15, 2
7, 17.5, 0
8, 17.5, 2
9, 20, 0
10, 20, 2
*ELEMENT, TYPE=CAX4, ELSET=S
1, 1, 3, 4, 2
2, 3, 5, 6, 4
3, 5, 7, 8, 6
4, 7, 9, 10, 8
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=S, MATERIAL=STEEL
*BOUNDARY
1, 2, 2
2, 2, 2
3, 2, 2
4, 2, 2
5, 2, 2
6, 2, 2
7, 2, 2
8, 2, 2
9, 2, 2
10, 2, 2
*STEP
*STATIC
*DLOAD
1, P4, 10
*NODE FILE
U
*EL FILE
S
*END STEP
`,
  },
  {
    id: "m3d4",
    name: "Membran M3D4",
    blurb: "Plane-Stress in 3D, ux=0.01",
    inp: `*HEADING
M3D4 membrane patch
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
*ELEMENT, TYPE=M3D4, ELSET=M
1, 1, 2, 3, 4
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.0
*MEMBRANE SECTION, ELSET=M, MATERIAL=STEEL
1.0
*BOUNDARY
1, 1, 1
4, 1, 1
1, 2, 2
2, 2, 2
1, 3, 3
2, 3, 3
3, 3, 3
4, 3, 3
2, 1, 1, 0.01
3, 1, 1, 0.01
*STEP
*STATIC
*NODE FILE
U
*EL FILE
S
*END STEP
`,
  },
  {
    id: "bar",
    name: "Zugstab fein",
    blurb: "6×2×2 Hexaeder, σ = 210",
    inp: tensionFine(),
  },
  {
    id: "beam2d",
    name: "Kragträger 2D",
    blurb: "CPS4, Vergleich Euler-Balken",
    inp: cantilever2d(),
  },
  {
    id: "cps8",
    name: "Kragträger CPS8",
    blurb: "Quadratische Vierecke, 2. Ordnung",
    inp: cantileverCps8(),
  },
  {
    id: "c3d20",
    name: "Zugstab C3D20",
    blurb: "Ein 20-Knoten-Hexaeder",
    inp: tensionC3d20.trimStart(),
  },
  {
    id: "beam3d",
    name: "Kragträger 3D",
    blurb: "C3D8-Balken, Last −Z",
    inp: cantilever3d(),
  },
  {
    id: "tets",
    name: "Tetraeder-Quader",
    blurb: "C3D4-Netz unter Zug",
    inp: tetBlock(),
  },
  {
    id: "hole",
    name: "Lochplatte",
    blurb: "Viertelmodell, Kerbspannung",
    inp: plateWithHole(),
  },
  {
    id: "nlgeom",
    name: "Fachwerk NLGEOM",
    blurb: "T3D2 *STEP, NLGEOM, u = FL/EA",
    inp: `*HEADING
T3D2 NLGEOM — kleine Dehnung, ux = FL/EA = 0.5
*NODE
1, 0, 0, 0
2, 1000, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
200
*BOUNDARY
1, 1, 3
2, 2, 3
*STEP, NLGEOM
*STATIC
*CLOAD
2, 1, 21000
*NODE FILE
U, RF
*EL FILE
S
*END STEP
`,
  },
  {
    id: "nlgeom-c3d8",
    name: "Hexaeder NLGEOM",
    blurb: "C3D8 Total-Lagrange, ux = 0.01",
    inp: `*HEADING
C3D8 NLGEOM — kleine Dehnung, ux = FL/EA = 0.01
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
5, 0, 0, 10
6, 10, 0, 10
7, 10, 10, 10
8, 0, 10, 10
*ELEMENT, TYPE=C3D8, ELSET=SOLID
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=SOLID, MATERIAL=STEEL
*NSET, NSET=FIXED
1, 4, 5, 8
*BOUNDARY
FIXED, 1, 1
1, 2, 3
4, 3, 3
*STEP, NLGEOM
*STATIC
*CLOAD
2, 1, 5250
3, 1, 5250
6, 1, 5250
7, 1, 5250
*NODE FILE
U, RF
*EL FILE
S
*END STEP
`,
  },
  {
    id: "nlgeom-stretch",
    name: "SVK-Zug λ=1.2",
    blurb: "C3D8 finite Dehnung, σxx=26.4",
    inp: `*HEADING
C3D8 St. Venant–Kirchhoff, λ=1.2, nu=0 → Cauchy σxx=26.4
*NODE
1, 0, 0, 0
2, 1, 0, 0
3, 1, 1, 0
4, 0, 1, 0
5, 0, 0, 1
6, 1, 0, 1
7, 1, 1, 1
8, 0, 1, 1
*ELEMENT, TYPE=C3D8, ELSET=S
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=STEEL
*ELASTIC
100, 0.0
*SOLID SECTION, ELSET=S, MATERIAL=STEEL
*BOUNDARY
1, 1, 3
4, 1, 1
4, 3, 3
5, 1, 2
8, 1, 1
2, 2, 3
3, 3, 3
6, 2, 2
2, 1, 1, 0.2
3, 1, 1, 0.2
6, 1, 1, 0.2
7, 1, 1, 0.2
*STEP, NLGEOM
*STATIC
*NODE FILE
U, RF
*EL FILE
S
*END STEP
`,
  },
  {
    id: "plastic",
    name: "Zugstab plastisch",
    blurb: "T3D2 *PLASTIC, σ = 250, u ≈ 3.1",
    inp: `*HEADING
T3D2 *PLASTIC — J2 1D, σ=250, ux ≈ 3.095
*NODE
1, 0, 0, 0
2, 1000, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*PLASTIC
210, 0.0
420, 0.01
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
200
*BOUNDARY
1, 1, 3
2, 2, 3
*STEP
*STATIC
*CLOAD
2, 1, 50000
*NODE FILE
U, RF
*EL FILE
S
*END STEP
`,
  },
  {
    id: "plastic-c3d8",
    name: "Hexaeder plastisch",
    blurb: "C3D8 J2, σ=250, ux≈0.031",
    inp: `*HEADING
C3D8 *PLASTIC J2, σ=250, ux ≈ 0.03095
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
5, 0, 0, 10
6, 10, 0, 10
7, 10, 10, 10
8, 0, 10, 10
*ELEMENT, TYPE=C3D8, ELSET=S
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*PLASTIC
210, 0.0
420, 0.01
*SOLID SECTION, ELSET=S, MATERIAL=STEEL
*NSET, NSET=FIXED
1, 4, 5, 8
*BOUNDARY
FIXED, 1, 1
1, 2, 3
4, 3, 3
*STEP
*STATIC
*CLOAD
2, 1, 6250
3, 1, 6250
6, 1, 6250
7, 1, 6250
*NODE FILE
U, RF
*EL FILE
S, PEEQ
*END STEP
`,
  },
  {
    id: "contact",
    name: "Kontakt zwei Würfel",
    blurb: "*CONTACT PAIR, Interface uz = −0.005",
    inp: `*HEADING
Zwei C3D8 in Serie über *CONTACT PAIR, uz_top = -0.01
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
5, 0, 0, 10
6, 10, 0, 10
7, 10, 10, 10
8, 0, 10, 10
9, 0, 0, 10
10, 10, 0, 10
11, 10, 10, 10
12, 0, 10, 10
13, 0, 0, 20
14, 10, 0, 20
15, 10, 10, 20
16, 0, 10, 20
*ELEMENT, TYPE=C3D8, ELSET=E1
1, 1, 2, 3, 4, 5, 6, 7, 8
*ELEMENT, TYPE=C3D8, ELSET=E2
2, 9, 10, 11, 12, 13, 14, 15, 16
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.0
*SOLID SECTION, ELSET=E1, MATERIAL=STEEL
*SOLID SECTION, ELSET=E2, MATERIAL=STEEL
*NSET, NSET=BOT
1, 2, 3, 4
*NSET, NSET=TOP
13, 14, 15, 16
*NSET, NSET=ALL
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16
*SURFACE, NAME=MASTER, TYPE=ELEMENT
1, S2
*SURFACE, NAME=SLAVE, TYPE=ELEMENT
2, S1
*SURFACE INTERACTION, NAME=INT
*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR
1e8
*CONTACT PAIR, INTERACTION=INT, TYPE=NODE TO SURFACE
SLAVE, MASTER
*BOUNDARY
ALL, 1, 2
BOT, 3, 3
TOP, 3, 3, -0.01
*STEP
*STATIC
*NODE FILE
U, RF
*EL FILE
S
*END STEP
`,
  },
  {
    id: "contact-friction",
    name: "Kontakt mit Reibung",
    blurb: "Coulomb μ=0.8, Interface haftet",
    inp: `*HEADING
Coulomb-Haften: Würfel auf Fundament, μ=0.8, ux_top=0.01
*NODE
1, 0, 0, -10
2, 10, 0, -10
3, 10, 10, -10
4, 0, 10, -10
5, 0, 0, 0
6, 10, 0, 0
7, 10, 10, 0
8, 0, 10, 0
9, 0, 0, -0.001
10, 10, 0, -0.001
11, 10, 10, -0.001
12, 0, 10, -0.001
13, 0, 0, 10
14, 10, 0, 10
15, 10, 10, 10
16, 0, 10, 10
*ELEMENT, TYPE=C3D8, ELSET=FND
1, 1, 2, 3, 4, 5, 6, 7, 8
*ELEMENT, TYPE=C3D8, ELSET=BLK
2, 9, 10, 11, 12, 13, 14, 15, 16
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.0
*SOLID SECTION, ELSET=FND, MATERIAL=STEEL
*SOLID SECTION, ELSET=BLK, MATERIAL=STEEL
*NSET, NSET=FOUND
1, 2, 3, 4, 5, 6, 7, 8
*NSET, NSET=TOP
13, 14, 15, 16
*NSET, NSET=ALL
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16
*SURFACE, NAME=MASTER, TYPE=ELEMENT
1, S2
*SURFACE, NAME=SLAVE, TYPE=ELEMENT
2, S1
*SURFACE INTERACTION, NAME=INT
*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR
1e8
*FRICTION
0.8
*CONTACT PAIR, INTERACTION=INT
SLAVE, MASTER
*BOUNDARY
FOUND, 1, 3
ALL, 2, 2
TOP, 1, 1, 0.01
TOP, 3, 3, -0.01
*STEP
*STATIC
*NODE FILE
U, RF
*EL FILE
S
*END STEP
`,
  },
  {
    id: "nlgeom-plastic",
    name: "NLGEOM + Plastizität",
    blurb: "C3D8 J2 mit NLGEOM, ux≈0.031",
    inp: `*HEADING
C3D8 *STEP, NLGEOM + *PLASTIC, σ=250, ux ≈ 0.03095
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 10, 10, 0
4, 0, 10, 0
5, 0, 0, 10
6, 10, 0, 10
7, 10, 10, 10
8, 0, 10, 10
*ELEMENT, TYPE=C3D8, ELSET=S
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*PLASTIC
210, 0.0
420, 0.01
*SOLID SECTION, ELSET=S, MATERIAL=STEEL
*NSET, NSET=FIXED
1, 4, 5, 8
*BOUNDARY
FIXED, 1, 1
1, 2, 3
4, 3, 3
*STEP, NLGEOM
*STATIC
*CLOAD
2, 1, 6250
3, 1, 6250
6, 1, 6250
7, 1, 6250
*NODE FILE
U, RF
*EL FILE
S, PEEQ
*END STEP
`,
  },
  {
    id: "heat",
    name: "Wärmeleitung T3D2",
    blurb: "*HEAT TRANSFER, T(L/2)=50",
    inp: `*HEADING
1D conduction, T(0)=0, T(L)=100 → Tmid=50
*NODE
1, 0, 0, 0
2, 500, 0, 0
3, 1000, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
2, 2, 3
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*CONDUCTIVITY
50
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
1
*BOUNDARY
1, 11, 11, 0
3, 11, 11, 100
*STEP
*HEAT TRANSFER, STEADY STATE
*NODE FILE
NT
*END STEP
`,
  },
  {
    id: "dynamic",
    name: "SDOF Newmark",
    blurb: "*DYNAMIC, freie Schwingung u(T/2)=−u0",
    inp: `*HEADING
T3D2 SDOF, ω=10 rad/s, u(0)=0.01 → u(π/10)=-0.01
*NODE
1, 0, 0, 0
2, 1, 0, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 2
*MATERIAL, NAME=STEEL
*ELASTIC
100, 0.0
*DENSITY
2.0
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
1.0
*BOUNDARY
1, 1, 3
2, 2, 3
*INITIAL CONDITIONS, TYPE=DISPLACEMENT
2, 1, 0.01
*STEP
*DYNAMIC
0.005, 0.3141592653589793
*NODE FILE
U
*END STEP
`,
  },
];

export const DEFAULT_EXAMPLE = EXAMPLES[0];
