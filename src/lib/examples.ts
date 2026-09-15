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
    id: "cax3",
    name: "Rohr CAX3",
    blurb: "Achsensymmetrie Dreieck, Innendruck",
    inp: `*HEADING
CAX3 thick cylinder
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
*ELEMENT, TYPE=CAX3, ELSET=S
1, 1, 3, 2
2, 3, 4, 2
3, 3, 5, 4
4, 5, 6, 4
5, 5, 7, 6
6, 7, 8, 6
7, 7, 9, 8
8, 9, 10, 8
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
1, P3, 10
*NODE FILE
U
*EL FILE
S
*END STEP
`,
  },
  {
    id: "c3d10t",
    name: "C3D10T B-bar",
    blurb: "Tetraeder mit mittlerer Dilatation, εx=0.001",
    inp: `*HEADING
C3D10T B-bar tet, confined εx=0.001
*NODE
1, 0, 0, 0
2, 10, 0, 0
3, 0, 10, 0
4, 0, 0, 10
5, 5, 0, 0
6, 5, 5, 0
7, 0, 5, 0
8, 0, 0, 5
9, 5, 0, 5
10, 0, 5, 5
*ELEMENT, TYPE=C3D10T, ELSET=S
1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=S, MATERIAL=STEEL
*BOUNDARY
1, 1, 3
2, 1, 1, 0.01
2, 2, 3
3, 1, 3
4, 1, 3
5, 1, 1, 0.005
5, 2, 3
6, 1, 1, 0.005
6, 2, 3
7, 1, 3
8, 1, 3
9, 1, 1, 0.005
9, 2, 3
10, 1, 3
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
    id: "bolted-joint",
    name: "Schraubenverbindung 30 kN",
    blurb: "C3D20-Lasche, Pretension 30 kN, Kontakt + TIE",
    inp: `** Generated by Mecway 33
*NODE
1,0.0375,0.005,0
2,0.05284746489505,0.01,0.025
3,0.057,0.01,0.025
4,0.0375,0.01,0
5,0.0375,0.005,0.05
6,0.0575806504495,0.01,0.02254032522475
7,0.05430545635174,0.01,0.01680545635174
8,0.0375,0.01,0.05
9,0.05861091270347,0.01,0.02111091270347
10,0.057,0.0075,0.025
11,0.05861091270347,0.0075,0.02111091270347
12,0.05284746489505,0.005,0.025
13,0.08125,0.01,0.0125
14,0.08141745418902,0.01,0.025
15,0.08125,0.015,0
16,0.08125,0.015,0.0125
17,0.057,0.005,0.025
18,0.0575806504495,0.005,0.02254032522475
19,0.05430545635174,0.005,0.01680545635174
20,0.05861091270347,0.005,0.02111091270347
21,0.0625,0.01,0.01534746489505
22,0.0375,0.005,0.025
23,0.0625,0.01,0.0195
24,0.0625,0.01,0.03465253510495
25,0.0625,0.01,0.0305
26,0.08141745418902,0.015,0.025
27,0.0625,0.0075,0.0195
28,0.0375,0.01,0.025
29,0.0625,0.005,0.01534746489505
30,0.0625,0.005,0.0195
31,0.0625,0.0075,0.0305
32,0.06638908729653,0.0075,0.02111091270347
33,0.0625,0.005,0.01117444782879
34,0.0625,0.005,0.0305
35,0.0625,0.0075,0
36,0.0625,0.0075,0.01117444909293
37,0.06004032522475,0.01,0.0299193495505
38,0.0625,0.01,0.005597464895051
39,0.0625,0.01,0.01117444782879
40,0.06004032522475,0.005,0.0299193495505
41,0.0625,0.005,0
42,0.0875,0.0075,0.0125
43,0.05430545635174,0.01,0.03319454364826
44,0.0875,0.0075,0.025
45,0.0875,0.005,0.0375
46,0.0625,0.01,0
47,0.0875,0.01,0.00625
48,0.05861091270347,0.01,0.02888908729653
49,0.0875,0.01,0.0125
50,0.0875,0.005,0.04375
51,0.0875,0.01,0.01875
52,0.0875,0.01,0.00625
53,0.0875,0.005,0
54,0.0875,0.01,0.0125
55,0.0875,0.005,0.025
56,0.0875,0.01,0.01875
57,0.0875,0.0125,0
58,0.0875,0.01,0
59,0.06638908729653,0.01,0.02111091270347
60,0.0875,0.01,0.025
61,0.07069454364826,0.01,0.03319454364826
62,0.06638908729653,0.01,0.02888908729653
63,0.0674193495505,0.01,0.02745967477525
64,0.06336824088834,0,0.02992403876506
65,0.05757596123494,0,0.02586824088833
66,0.0625,0.005,0.03465253510495
67,0.06495967477525,0.01,0.0200806504495
68,0.05877175316372,0.015,0.03032448828788
69,0.05791138850919,0.015,0.03155321635431
70,0.06782448828788,0.015,0.02872824683628
71,0.06905321635431,0.015,0.02958861149081
72,0.06782448828788,0.02,0.02872824683628
73,0.06905321635431,0.02,0.02958861149081
74,0.05877175316372,0.02,0.03032448828788
75,0.05791138850919,0.02,0.03155321635431
76,0.05717551171212,0.015,0.02127175316372
77,0.05594678364569,0.015,0.02041138850919
78,0.05717551171212,0.02,0.02127175316372
79,0.05594678364569,0.02,0.02041138850919
80,0.06622824683628,0.015,0.01967551171212
81,0.06708861149081,0.015,0.01844678364569
82,0.06622824683628,0.02,0.01967551171212
83,0.06708861149081,0.02,0.01844678364569
84,0.05717551171212,0.005,0.02127175316372
85,0.05594678364569,0.005,0.02041138850919
86,0.06622824683628,0.005,0.01967551171212
87,0.06708861149081,0.005,0.01844678364569
88,0.06622824683628,0,0.01967551171212
89,0.06708861149081,0,0.01844678364569
90,0.05717551171212,0,0.02127175316372
91,0.05594678364569,0,0.02041138850919
92,0.05877175316372,0.005,0.03032448828788
93,0.05791138850919,0.005,0.03155321635431
94,0.05877175316372,0,0.03032448828788
95,0.05791138850919,0,0.03155321635431
96,0.06782448828788,0.005,0.02872824683628
97,0.09375,0.01,0.0125
98,0.1,0.01,0.0125
99,0.06875,0.005,0.03815253510495
100,0.10625,0.01,0.0125
101,0.06875,0.005,0.05
102,0.1125,0.01,0.0125
103,0.06875,0.01,0.03815253510495
104,0.06875,0.01,0.05
105,0.11875,0.01,0.0125
106,0.07565253510495,0.005,0.03125
107,0.075,0.005,0.0375
108,0.075,0.005,0.04375
109,0.075,0.005,0.05
110,0.125,0.01,0.0125
111,0.09375,0.015,0.0125
112,0.075,0.0075,0.0375
113,0.075,0.0075,0.05
114,0.1,0.015,0.0125
115,0.07565253510495,0.01,0.03125
116,0.075,0.01,0.0375
117,0.075,0.01,0.04375
118,0.075,0.01,0.05
119,0.10625,0.015,0.0125
120,0.08125,0.005,0.0375
121,0.0625,0.0075,0.03882555090707
122,0.08125,0.005,0.05
123,0.0625,0.0075,0.05
124,0.0575806504495,0.01,0.02745967477525
125,0.1125,0.015,0.0125
126,0.0625,0.01,0.03882555217121
127,0.08125,0.01,0.0375
128,0.0625,0.01,0.04440253510495
129,0.08125,0.01,0.05
130,0.0875,0.0075,0.0375
131,0.0625,0.005,0.05
132,0.0875,0.005,0.03125
133,0.0875,0.0075,0.05
134,0.11875,0.015,0.0125
135,0.0875,0.01,0.03125
136,0.0625,0.01,0.05
137,0.125,0.015,0.0125
138,0.0875,0.01,0.0375
139,0.09375,0.015,0.025
140,0.0875,0.01,0.04375
141,0.1,0.015,0.025
142,0.10625,0.015,0.025
143,0.0875,0.005,0.05
144,0.1125,0.015,0.025
145,0.11875,0.015,0.025
146,0.125,0.015,0.025
147,0.09375,0.01,0.025
148,0.0875,0.01,0.05
149,0.1,0.01,0.025
150,0.10625,0.01,0.025
151,0.1125,0.01,0.025
152,0.11875,0.01,0.025
153,0.125,0.01,0.025
154,0.09375,0.01,0
155,0.1,0.01,0
156,0.10625,0.01,0
157,0.1125,0.01,0
158,0.11875,0.01,0
159,0.125,0.01,0
160,0.09375,0.015,0
161,0.1,0.015,0
162,0.07215253510495,0.01,0.025
163,0.068,0.01,0.025
164,0.0674193495505,0.01,0.02254032522475
165,0.07069454364826,0.01,0.01680545635174
166,0.0375,0.005,0.03125
167,0.0375,0.005,0.0375
168,0.0375,0.005,0.04375
169,0.0375,0.0075,0.025
170,0.0375,0.0075,0.0375
171,0.0375,0.0075,0.05
172,0.0375,0.01,0.03125
173,0.0375,0.01,0.0375
174,0.0375,0.01,0.04375
175,0.04309746489505,0.005,0.025
176,0.04375,0.005,0.0375
177,0.04375,0.005,0.05
178,0.04309746489505,0.01,0.025
179,0.04375,0.01,0.0375
180,0.04375,0.01,0.05
181,0.06638908729653,0.01,0.02111091270347
182,0.04934746489505,0.005,0.03125
183,0.05,0.005,0.0375
184,0.05,0.005,0.04375
185,0.05,0.005,0.05
186,0.10625,0.015,0
187,0.05,0.0075,0.0375
188,0.05,0.0075,0.05
189,0.1125,0.015,0
190,0.04934746489505,0.01,0.03125
191,0.05,0.01,0.0375
192,0.05,0.01,0.04375
193,0.05,0.01,0.05
194,0.11875,0.015,0
195,0.05625,0.005,0.03815253510495
196,0.05625,0.005,0.05
197,0.125,0.015,0
198,0.05625,0.01,0.03815253510495
199,0.05625,0.01,0.05
200,0.05861091270347,0.0075,0.02888908729653
201,0.0625,0.005,0.03882555217121
202,0.0625,0.005,0.04440253510495
203,0.068,0.0075,0.025
204,0.0375,0.005,0.00625
205,0.0375,0.005,0.0125
206,0.0375,0.005,0.01875
207,0.0375,0.0075,0
208,0.0375,0.0075,0.0125
209,0.09375,0.015,0.0375
210,0.0375,0.01,0.00625
211,0.0375,0.01,0.0125
212,0.0375,0.01,0.01875
213,0.04375,0.005,0
214,0.04375,0.005,0.0125
215,0.1,0.015,0.0375
216,0.04375,0.01,0
217,0.04375,0.01,0.0125
218,0.07069454364826,0.01,0.03319454364826
219,0.05,0.005,0
220,0.05,0.005,0.00625
221,0.05,0.005,0.0125
222,0.04934746489505,0.005,0.01875
223,0.04867444782879,0.005,0.025
224,0.05,0.0075,0
225,0.05,0.0075,0.0125
226,0.04867444909293,0.0075,0.025
227,0.05,0.01,0
228,0.05,0.01,0.00625
229,0.05,0.01,0.0125
230,0.04934746489505,0.01,0.01875
231,0.04867444782879,0.01,0.025
232,0.05625,0.005,0
233,0.05625,0.005,0.01184746489505
234,0.06495967477525,0.01,0.0200806504495
235,0.05625,0.01,0
236,0.05625,0.01,0.01184746489505
237,0.05430545635174,0.005,0.03319454364826
238,0.0625,0.005,0.005597464895051
239,0.10625,0.015,0.0375
240,0.1125,0.015,0.0375
241,0.11875,0.015,0.0375
242,0.06638908729653,0.01,0.02888908729653
243,0.0674193495505,0.01,0.02745967477525
244,0.07215253510495,0.005,0.025
245,0.068,0.005,0.025
246,0.0674193495505,0.005,0.02254032522475
247,0.07069454364826,0.005,0.01680545635174
248,0.06875,0.005,0
249,0.06875,0.005,0.01184746489505
250,0.05861091270347,0.005,0.02888908729653
251,0.06875,0.01,0
252,0.06875,0.01,0.01184746489505
253,0.0575806504495,0.005,0.02745967477525
254,0.075,0.005,0
255,0.075,0.005,0.00625
256,0.075,0.005,0.0125
257,0.07565253510495,0.005,0.01875
258,0.07632555217121,0.005,0.025
259,0.075,0.0075,0
260,0.075,0.0075,0.0125
261,0.07632555092178,0.0075,0.02499999986443
262,0.075,0.01,0
263,0.075,0.01,0.00625
264,0.075,0.01,0.0125
265,0.07565253510495,0.01,0.01875
266,0.07632555217121,0.01,0.025
267,0.08125,0.005,0
268,0.08125,0.005,0.0125
269,0.08141745418902,0.005,0.025
270,0.08125,0.01,0
271,0.08125,0.01,0.0125
272,0.08141745418902,0.01,0.025
273,0.0875,0.005,0.00625
274,0.0875,0.005,0.0125
275,0.0875,0.005,0.01875
276,0.0875,0.0075,0
277,0.06638908729653,0.005,0.02111091270347
278,0.07069454364826,0.005,0.03319454364826
279,0.06638908729653,0.005,0.02888908729653
280,0.0674193495505,0.005,0.02745967477525
281,0.125,0.015,0.0375
282,0.09375,0.01,0.0375
283,0.06495967477525,0.005,0.0200806504495
284,0.1,0.01,0.0375
285,0.10625,0.01,0.0375
286,0.1125,0.01,0.0375
287,0.11875,0.01,0.0375
288,0.125,0.01,0.0375
289,0.09375,0.015,0.05
290,0.1,0.015,0.05
291,0.10625,0.015,0.05
292,0.1125,0.015,0.05
293,0.11875,0.015,0.05
294,0.125,0.015,0.05
295,0.09375,0.01,0.05
296,0.1,0.01,0.05
297,0.10625,0.01,0.05
298,0.1125,0.01,0.05
299,0.11875,0.01,0.05
300,0.125,0.01,0.05
301,0.1,0.0125,0.0125
302,0.1125,0.0125,0.0125
303,0.125,0.0125,0.0125
304,0.1,0.015,0.01875
305,0.1125,0.015,0.01875
306,0.125,0.015,0.01875
307,0.1,0.0125,0.025
308,0.1125,0.0125,0.025
309,0.125,0.0125,0.025
310,0.1,0.01,0.01875
311,0.1125,0.01,0.01875
312,0.125,0.01,0.01875
313,0.1,0.0125,0
314,0.1125,0.0125,0
315,0.06495967477525,0.01,0.0299193495505
316,0.125,0.0125,0
317,0.1,0.015,0.00625
318,0.1125,0.015,0.00625
319,0.125,0.015,0.00625
320,0.06495967477525,0.005,0.0299193495505
321,0.1,0.01,0.00625
322,0.1125,0.01,0.00625
323,0.06638908729653,0.0075,0.02888908729653
324,0.125,0.01,0.00625
325,0.1,0.015,0.03125
326,0.1125,0.015,0.03125
327,0.125,0.015,0.03125
328,0.1,0.0125,0.0375
329,0.1125,0.0125,0.0375
330,0.125,0.0125,0.0375
331,0.1,0.01,0.03125
332,0.1125,0.01,0.03125
333,0.125,0.01,0.03125
334,0.1,0.015,0.04375
335,0.1125,0.015,0.04375
336,0.125,0.015,0.04375
337,0.1,0.0125,0.05
338,0.1125,0.0125,0.05
339,0.125,0.0125,0.05
340,0.1,0.01,0.04375
341,0.1125,0.01,0.04375
342,0.125,0.01,0.04375
343,0.06004032522475,0.01,0.0200806504495
344,0.06004032522475,0.005,0.0200806504495
345,0.0375,0.01,0
346,0.05284746489505,0.015,0.025
347,0.057,0.015,0.025
348,0.0375,0.015,0
349,0.0375,0.01,0.05
350,0.0575806504495,0.015,0.02254032522475
351,0.05430545635174,0.015,0.01680545635174
352,0.0375,0.015,0.05
353,0.05861091270347,0.015,0.02111091270347
354,0.057,0.0125,0.025
355,0.05861091270347,0.0125,0.02111091270347
356,0.05284746489505,0.01,0.025
357,0.06905321635431,0.005,0.02958861149081
358,0.06782448828788,0,0.02872824683628
359,0.06495967477525,0.015,0.0299193495505
360,0.06905321635431,0,0.02958861149081
361,0.057,0.01,0.025
362,0.0575806504495,0.01,0.02254032522475
363,0.05430545635174,0.01,0.01680545635174
364,0.05861091270347,0.01,0.02111091270347
365,0.0625,0.015,0.01534746489505
366,0.0375,0.01,0.025
367,0.0625,0.015,0.0195
368,0.0625,0.015,0.03465253510495
369,0.0625,0.015,0.0305
370,0.06388918542134,0.015,0.0328784620241
371,0.0625,0.0125,0.0195
372,0.0375,0.015,0.025
373,0.0625,0.01,0.01534746489505
374,0.0625,0.01,0.0195
375,0.0625,0.0125,0.0305
376,0.06638908729653,0.0125,0.02111091270347
377,0.0625,0.01,0.01117444782879
378,0.0625,0.01,0.0305
379,0.0625,0.0125,0
380,0.0625,0.0125,0.01117444909293
381,0.06004032522475,0.015,0.0299193495505
382,0.0625,0.015,0.005597464895051
383,0.0625,0.015,0.01117444782879
384,0.06004032522475,0.01,0.0299193495505
385,0.0625,0.01,0
386,0.0875,0.0125,0.0125
387,0.05430545635174,0.015,0.03319454364826
388,0.0875,0.0125,0.025
389,0.0875,0.01,0.0375
390,0.0625,0.015,0
391,0.0875,0.015,0.00625
392,0.05861091270347,0.015,0.02888908729653
393,0.0875,0.015,0.0125
394,0.0875,0.01,0.04375
395,0.0875,0.015,0.01875
396,0.06905321635431,0.0175,0.02958861149081
397,0.0875,0.01,0
398,0.06388918542134,0.02,0.0328784620241
399,0.0875,0.01,0.025
400,0.06495967477525,0.01,0.0299193495505
401,0.05791138850919,0.0175,0.03155321635431
402,0.0875,0.015,0
403,0.0546215379759,0.015,0.02638918542134
404,0.0875,0.015,0.025
405,0.06638908729653,0.0125,0.02888908729653
406,0.0546215379759,0.02,0.02638918542134
407,0.05594678364569,0.0175,0.02041138850919
408,0.0703784620241,0.015,0.02361081457866
409,0.06708861149081,0.0175,0.01844678364569
410,0.0625,0.01,0.03465253510495
411,0.0703784620241,0.02,0.02361081457866
412,0.06111081457867,0.015,0.0171215379759
413,0.06111081457867,0.02,0.0171215379759
414,0.06111081457867,0.005,0.0171215379759
415,0.06708861149081,0.0025,0.01844678364569
416,0.06111081457867,0,0.0171215379759
417,0.05594678364569,0.0025,0.02041138850919
418,0.0546215379759,0.005,0.02638918542134
419,0.0546215379759,0,0.02638918542134
420,0.05791138850919,0.0025,0.03155321635431
421,0.06388918542134,0.005,0.0328784620241
422,0.06388918542134,0,0.0328784620241
423,0.06905321635431,0.0025,0.02958861149081
424,0.0703784620241,0.005,0.02361081457866
425,0.0703784620241,0,0.02361081457866
426,0.06004032522475,0.015,0.0200806504495
427,0.06004032522475,0.01,0.0200806504495
428,0.05840423977856,0.005,0.02213211781824
429,0.06536788218176,0.005,0.02090423977856
430,0.06659576022145,0.005,0.02786788218176
431,0.05963211781824,0.005,0.02909576022145
432,0.06163175911167,0.005,0.02007596123494
433,0.06742403876506,0.005,0.02413175911167
434,0.06336824088834,0.005,0.02992403876506
435,0.05757596123494,0.005,0.02586824088833
436,0.05840423977856,0.0075,0.02213211781824
437,0.05840423977856,0.01,0.02213211781824
438,0.05840423977856,0.0125,0.02213211781824
439,0.05840423977856,0.015,0.02213211781824
440,0.05963211781824,0.0075,0.02909576022145
441,0.05963211781824,0.01,0.02909576022145
442,0.05963211781824,0.0125,0.02909576022145
443,0.06875,0.01,0.03815253510495
444,0.05963211781824,0.015,0.02909576022145
445,0.06875,0.01,0.05
446,0.06659576022145,0.0075,0.02786788218176
447,0.06875,0.015,0.03815253510495
448,0.06875,0.015,0.05
449,0.06659576022145,0.01,0.02786788218176
450,0.07565253510495,0.01,0.03125
451,0.075,0.01,0.0375
452,0.075,0.01,0.04375
453,0.075,0.01,0.05
454,0.06659576022145,0.0125,0.02786788218176
455,0.06659576022145,0.015,0.02786788218176
456,0.075,0.0125,0.0375
457,0.075,0.0125,0.05
458,0.06536788218176,0.0075,0.02090423977856
459,0.07565253510495,0.015,0.03125
460,0.075,0.015,0.0375
461,0.075,0.015,0.04375
462,0.075,0.015,0.05
463,0.06536788218176,0.01,0.02090423977856
464,0.08125,0.01,0.0375
465,0.0625,0.0125,0.03882555090707
466,0.08125,0.01,0.05
467,0.0625,0.0125,0.05
468,0.0575806504495,0.015,0.02745967477525
469,0.06536788218176,0.0125,0.02090423977856
470,0.0625,0.015,0.03882555217121
471,0.08125,0.015,0.0375
472,0.0625,0.015,0.04440253510495
473,0.08125,0.015,0.05
474,0.0875,0.0125,0.0375
475,0.0625,0.01,0.05
476,0.0875,0.01,0.03125
477,0.0875,0.0125,0.05
478,0.06536788218176,0.015,0.02090423977856
479,0.0875,0.015,0.03125
480,0.0625,0.015,0.05
481,0.05757596123494,0.01,0.02586824088833
482,0.0875,0.015,0.0375
483,0.05757596123494,0.015,0.02586824088833
484,0.0875,0.015,0.04375
485,0.06336824088834,0.01,0.02992403876506
486,0.06336824088834,0.015,0.02992403876506
487,0.0875,0.01,0.05
488,0.06742403876506,0.01,0.02413175911167
489,0.06742403876506,0.015,0.02413175911167
490,0.06163175911167,0.01,0.02007596123494
491,0.06163175911167,0.015,0.02007596123494
492,0.0875,0.015,0.05
493,0.05840423977856,0.0175,0.02213211781824
494,0.05840423977856,0.02,0.02213211781824
495,0.05963211781824,0.0175,0.02909576022145
496,0.05963211781824,0.02,0.02909576022145
497,0.06659576022145,0.0175,0.02786788218176
498,0.06659576022145,0.02,0.02786788218176
499,0.06536788218176,0.0175,0.02090423977856
500,0.06536788218176,0.02,0.02090423977856
501,0.05757596123494,0.02,0.02586824088833
502,0.06336824088834,0.02,0.02992403876506
503,0.06742403876506,0.02,0.02413175911167
504,0.06163175911167,0.02,0.02007596123494
505,0.05840423977856,0.0025,0.02213211781824
506,0.07215253510495,0.015,0.025
507,0.068,0.015,0.025
508,0.0674193495505,0.015,0.02254032522475
509,0.07069454364826,0.015,0.01680545635174
510,0.0375,0.01,0.03125
511,0.0375,0.01,0.0375
512,0.0375,0.01,0.04375
513,0.0375,0.0125,0.025
514,0.0375,0.0125,0.0375
515,0.0375,0.0125,0.05
516,0.0375,0.015,0.03125
517,0.0375,0.015,0.0375
518,0.0375,0.015,0.04375
519,0.04309746489505,0.01,0.025
520,0.04375,0.01,0.0375
521,0.04375,0.01,0.05
522,0.04309746489505,0.015,0.025
523,0.04375,0.015,0.0375
524,0.04375,0.015,0.05
525,0.06638908729653,0.015,0.02111091270347
526,0.04934746489505,0.01,0.03125
527,0.05,0.01,0.0375
528,0.05,0.01,0.04375
529,0.05,0.01,0.05
530,0.05840423977856,0,0.02213211781824
531,0.05,0.0125,0.0375
532,0.05,0.0125,0.05
533,0.06536788218176,0.0025,0.02090423977856
534,0.04934746489505,0.015,0.03125
535,0.05,0.015,0.0375
536,0.05,0.015,0.04375
537,0.05,0.015,0.05
538,0.06536788218176,0,0.02090423977856
539,0.05625,0.01,0.03815253510495
540,0.05625,0.01,0.05
541,0.06659576022145,0.0025,0.02786788218176
542,0.05625,0.015,0.03815253510495
543,0.05625,0.015,0.05
544,0.05861091270347,0.0125,0.02888908729653
545,0.0625,0.01,0.03882555217121
546,0.0625,0.01,0.04440253510495
547,0.068,0.0125,0.025
548,0.0375,0.01,0.00625
549,0.0375,0.01,0.0125
550,0.0375,0.01,0.01875
551,0.0375,0.0125,0
552,0.0375,0.0125,0.0125
553,0.06659576022145,0,0.02786788218176
554,0.0375,0.015,0.00625
555,0.0375,0.015,0.0125
556,0.0375,0.015,0.01875
557,0.04375,0.01,0
558,0.04375,0.01,0.0125
559,0.05963211781824,0.0025,0.02909576022145
560,0.04375,0.015,0
561,0.04375,0.015,0.0125
562,0.07069454364826,0.015,0.03319454364826
563,0.05,0.01,0
564,0.05,0.01,0.00625
565,0.05,0.01,0.0125
566,0.04934746489505,0.01,0.01875
567,0.04867444782879,0.01,0.025
568,0.05,0.0125,0
569,0.05,0.0125,0.0125
570,0.04867444909293,0.0125,0.025
571,0.05,0.015,0
572,0.05,0.015,0.00625
573,0.05,0.015,0.0125
574,0.04934746489505,0.015,0.01875
575,0.04867444782879,0.015,0.025
576,0.05625,0.01,0
577,0.05625,0.01,0.01184746489505
578,0.06495967477525,0.015,0.0200806504495
579,0.05625,0.015,0
580,0.05625,0.015,0.01184746489505
581,0.05430545635174,0.01,0.03319454364826
582,0.0625,0.01,0.005597464895051
583,0.05963211781824,0,0.02909576022145
584,0.06163175911167,0,0.02007596123494
585,0.06742403876506,0,0.02413175911167
586,0.06638908729653,0.015,0.02888908729653
587,0.0674193495505,0.015,0.02745967477525
588,0.07215253510495,0.01,0.025
589,0.068,0.01,0.025
590,0.0674193495505,0.01,0.02254032522475
591,0.07069454364826,0.01,0.01680545635174
592,0.06875,0.01,0
593,0.06875,0.01,0.01184746489505
594,0.05861091270347,0.01,0.02888908729653
595,0.06875,0.015,0
596,0.06875,0.015,0.01184746489505
597,0.0575806504495,0.01,0.02745967477525
598,0.075,0.01,0
599,0.075,0.01,0.00625
600,0.075,0.01,0.0125
601,0.07565253510495,0.01,0.01875
602,0.07632555217121,0.01,0.025
603,0.075,0.0125,0
604,0.075,0.0125,0.0125
605,0.07632555092178,0.0125,0.02499999986443
606,0.075,0.015,0
607,0.075,0.015,0.00625
608,0.075,0.015,0.0125
609,0.07565253510495,0.015,0.01875
610,0.07632555217121,0.015,0.025
611,0.08125,0.01,0
612,0.03125,0.01,0.025
613,0.025,0.01,0.025
614,0.01875,0.01,0.025
615,0.0125,0.01,0.025
616,0.00625,0.01,0.025
617,0,0.01,0.025
618,0.03125,0.005,0.025
619,0.025,0.005,0.025
620,0.01875,0.005,0.025
621,0.0125,0.005,0.025
622,0.00625,0.005,0.025
623,0,0.005,0.025
624,0.03125,0.005,0.0375
625,0.025,0.005,0.0375
626,0.01875,0.005,0.0375
627,0.0125,0.005,0.0375
628,0.00625,0.005,0.0375
629,0,0.005,0.0375
630,0.03125,0.01,0.0375
631,0.025,0.01,0.0375
632,0.01875,0.01,0.0375
633,0.0125,0.01,0.0375
634,0.00625,0.01,0.0375
635,0,0.01,0.0375
636,0.03125,0.01,0.0125
637,0.025,0.01,0.0125
638,0.01875,0.01,0.0125
639,0.0125,0.01,0.0125
640,0.00625,0.01,0.0125
641,0,0.01,0.0125
642,0.03125,0.005,0.0125
643,0.025,0.005,0.0125
644,0.01875,0.005,0.0125
645,0.0125,0.005,0.0125
646,0.00625,0.005,0.0125
647,0,0.005,0.0125
648,0.03125,0.005,0.05
649,0.025,0.005,0.05
650,0.01875,0.005,0.05
651,0.0125,0.005,0.05
652,0.00625,0.005,0.05
653,0,0.005,0.05
654,0.03125,0.01,0.05
655,0.025,0.01,0.05
656,0.01875,0.01,0.05
657,0.0125,0.01,0.05
658,0.00625,0.01,0.05
659,0,0.01,0.05
660,0.03125,0.01,0
661,0.025,0.01,0
662,0.01875,0.01,0
663,0.0125,0.01,0
664,0.00625,0.01,0
665,0,0.01,0
666,0.03125,0.005,0
667,0.025,0.005,0
668,0.01875,0.005,0
669,0.0125,0.005,0
670,0.00625,0.005,0
671,0,0.005,0
672,0.025,0.0075,0.025
673,0.0125,0.0075,0.025
674,0,0.0075,0.025
675,0.025,0.005,0.03125
676,0.0125,0.005,0.03125
677,0,0.005,0.03125
678,0.025,0.0075,0.0375
679,0.0125,0.0075,0.0375
680,0,0.0075,0.0375
681,0.025,0.01,0.03125
682,0.0125,0.01,0.03125
683,0,0.01,0.03125
684,0.025,0.0075,0.0125
685,0.0125,0.0075,0.0125
686,0,0.0075,0.0125
687,0.025,0.005,0.01875
688,0.0125,0.005,0.01875
689,0,0.005,0.01875
690,0.025,0.01,0.01875
691,0.0125,0.01,0.01875
692,0,0.01,0.01875
693,0.025,0.005,0.04375
694,0.0125,0.005,0.04375
695,0,0.005,0.04375
696,0.025,0.0075,0.05
697,0.0125,0.0075,0.05
698,0,0.0075,0.05
699,0.025,0.01,0.04375
700,0.0125,0.01,0.04375
701,0,0.01,0.04375
702,0.025,0.0075,0
703,0.0125,0.0075,0
704,0,0.0075,0
705,0.025,0.005,0.00625
706,0.0125,0.005,0.00625
707,0,0.005,0.00625
708,0.025,0.01,0.00625
709,0.0125,0.01,0.00625
710,0,0.01,0.00625
711,0.05840423977856,0.009875,0.02213211781824
712,0.05840423977856,0.00975,0.02213211781824
713,0.05963211781824,0.009875,0.02909576022145
714,0.05963211781824,0.00975,0.02909576022145
715,0.06659576022145,0.009875,0.02786788218176
716,0.06659576022145,0.00975,0.02786788218176
717,0.06536788218176,0.009875,0.02090423977856
718,0.06536788218176,0.00975,0.02090423977856
719,0.05757596123494,0.00975,0.02586824088833
720,0.06336824088834,0.00975,0.02992403876506
721,0.06742403876506,0.00975,0.02413175911167
722,0.06163175911167,0.00975,0.02007596123494
723,0,0,0
*ELEMENT,TYPE=C3D20
1,221,33,39,229,20,30,23,9,233,36,236,225,344,27,343,
11,19,29,21,7
2,478,439,494,500,81,77,79,83,491,493,504,499,412,407,413,
409,80,76,78,82
3,33,256,264,39,30,277,181,23,249,260,252,36,283,32,234,
27,29,247,165,21
4,256,274,49,264,258,55,60,266,268,42,271,260,269,44,272,
261,257,275,51,265
5,428,429,538,530,85,87,89,91,432,533,584,505,414,415,416,
417,84,86,88,90
6,431,428,530,583,93,85,91,95,435,505,65,559,418,417,419,
420,92,84,90,94
7,527,535,470,545,594,392,369,378,531,542,465,539,544,381,375,
384,581,387,368,410
8,54,393,404,399,98,114,141,149,386,395,388,56,301,304,307,
310,97,111,139,147
9,201,107,116,126,131,109,118,136,99,112,103,121,101,113,104,
123,202,108,117,128
10,183,201,126,191,185,131,136,193,195,121,198,187,196,123,199,
188,184,202,128,192
11,107,45,138,116,109,143,148,118,120,130,127,112,122,133,129,
113,108,50,140,117
12,98,114,141,149,102,125,144,151,301,304,307,310,302,305,308,
311,100,119,142,150
13,266,258,107,116,163,245,279,242,261,106,112,115,203,280,323,
243,162,244,278,218
14,221,229,231,223,20,9,3,17,225,230,226,222,11,6,10,
18,19,7,2,12
15,102,125,144,151,110,137,146,153,302,305,308,311,303,306,309,
312,105,134,145,152
16,397,402,393,54,155,161,114,98,57,391,386,52,313,317,301,
321,154,160,111,97
17,22,223,231,28,167,183,191,173,175,226,178,169,176,187,179,
170,166,182,190,172
18,167,183,191,173,5,185,193,8,176,187,179,170,177,188,180,
171,168,184,192,174
19,264,256,258,266,181,277,245,163,260,257,261,265,32,246,203,
164,165,247,244,162
20,1,219,227,4,205,221,229,211,213,224,216,207,214,225,217,
208,204,220,228,210
21,205,221,229,211,22,223,231,28,214,225,217,208,175,226,178,
169,206,222,230,212
22,219,41,46,227,221,33,39,229,232,35,235,224,233,36,236,
225,220,238,38,228
23,41,254,262,46,33,256,264,39,248,259,251,35,249,260,252,
36,238,255,263,38
24,201,126,116,107,34,25,242,279,121,103,112,99,31,315,323,
320,66,24,218,278
25,254,53,58,262,256,274,49,264,267,276,270,259,268,42,271,
260,255,273,47,263
26,155,161,114,98,157,189,125,102,313,317,301,321,314,318,302,
322,156,186,119,100
27,157,189,125,102,159,197,137,110,314,318,302,322,316,319,303,
324,158,194,134,105
28,399,404,482,389,149,141,215,284,388,479,474,476,307,325,328,
331,147,139,209,282
29,149,141,215,284,151,144,240,286,307,325,328,331,308,326,329,
332,150,142,239,285
30,151,144,240,286,153,146,281,288,308,326,329,332,309,327,330,
333,152,145,241,287
31,389,482,492,487,284,215,290,296,474,484,477,394,328,334,337,
340,282,209,289,295
32,183,191,126,201,250,48,25,34,187,198,121,195,200,37,31,
40,237,43,24,66
33,223,231,191,183,17,3,48,250,226,190,187,182,10,124,200,
253,12,2,43,237
34,284,215,290,296,286,240,292,298,328,334,337,340,329,335,338,
341,285,239,291,297
35,258,55,60,266,107,45,138,116,269,44,272,261,120,130,127,
112,106,132,135,115
36,286,240,292,298,288,281,294,300,329,335,338,341,330,336,339,
342,287,241,293,299
37,565,377,383,573,364,374,367,353,577,380,580,569,427,371,426,
355,363,373,365,351
38,567,575,535,527,361,347,392,594,570,534,531,526,354,468,544,
597,356,346,387,581
39,377,600,608,383,374,59,525,367,593,604,596,380,67,376,578,
371,373,591,509,365
40,600,54,393,608,602,399,404,610,13,386,16,604,14,388,26,
605,601,56,395,609
41,430,431,583,553,357,93,95,360,434,559,64,541,421,420,422,
423,96,92,94,358
42,602,399,404,610,451,389,482,460,14,388,26,605,464,474,471,
456,450,476,479,459
43,429,430,553,538,87,357,360,89,433,541,585,533,424,423,425,
415,86,96,358,88
44,437,441,449,463,439,444,455,478,481,485,488,490,483,486,489,
491,438,442,454,469
45,545,451,460,470,475,453,462,480,443,456,447,465,445,457,448,
467,546,452,461,472
46,527,545,470,535,529,475,480,537,539,465,542,531,540,467,543,
532,528,546,472,536
47,451,389,482,460,453,487,492,462,464,474,471,456,466,477,473,
457,452,394,484,461
48,428,431,430,429,712,714,716,718,435,434,433,432,719,720,721,
722,436,440,446,458
49,610,602,451,460,507,589,62,586,605,450,456,459,547,63,405,
587,506,588,61,562
50,565,573,575,567,364,353,347,361,569,574,570,566,355,350,354,
362,363,351,346,356
51,439,444,455,478,494,496,498,500,483,486,489,491,501,502,503,
504,493,495,497,499
52,428,429,430,431,530,538,553,583,432,433,434,435,584,585,64,
65,505,533,541,559
53,366,567,575,372,511,527,535,517,519,570,522,513,520,531,523,
514,510,526,534,516
54,511,527,535,517,349,529,537,352,520,531,523,514,521,532,524,
515,512,528,536,518
55,608,600,602,610,525,59,589,507,604,601,605,609,376,590,547,
508,509,591,588,506
56,345,563,571,348,549,565,573,555,557,568,560,551,558,569,561,
552,548,564,572,554
57,549,565,573,555,366,567,575,372,558,569,561,552,519,570,522,
513,550,566,574,556
58,563,385,390,571,565,377,383,573,576,379,579,568,577,380,580,
569,564,582,382,572
59,385,598,606,390,377,600,608,383,592,603,595,379,593,604,596,
380,582,599,607,382
60,545,470,460,451,378,369,586,62,465,447,456,443,375,359,405,
400,410,368,562,61
61,598,397,402,606,600,54,393,608,611,57,15,603,13,386,16,
604,599,52,391,607
62,444,455,498,496,69,71,73,75,486,497,502,495,370,396,398,
401,68,70,72,74
63,439,444,496,494,77,69,75,79,483,495,501,493,403,401,406,
407,76,68,74,78
64,455,478,500,498,71,81,83,73,489,499,503,497,408,409,411,
396,70,80,82,72
65,28,22,167,173,613,619,625,631,169,166,170,172,672,675,678,
681,612,618,624,630
66,613,619,625,631,615,621,627,633,672,675,678,681,673,676,679,
682,614,620,626,632
67,615,621,627,633,617,623,629,635,673,676,679,682,674,677,680,
683,616,622,628,634
68,211,205,22,28,637,643,619,613,208,206,169,212,684,687,672,
690,636,642,618,612
69,637,643,619,613,639,645,621,615,684,687,672,690,685,688,673,
691,638,644,620,614
70,639,645,621,615,641,647,623,617,685,688,673,691,686,689,674,
692,640,646,622,616
71,173,167,5,8,631,625,649,655,170,168,171,174,678,693,696,
699,630,624,648,654
72,631,625,649,655,633,627,651,657,678,693,696,699,679,694,697,
700,632,626,650,656
73,633,627,651,657,635,629,653,659,679,694,697,700,680,695,698,
701,634,628,652,658
74,4,1,205,211,661,667,643,637,207,204,208,210,702,705,684,
708,660,666,642,636
75,661,667,643,637,663,669,645,639,702,705,684,708,703,706,685,
709,662,668,644,638
76,663,669,645,639,665,671,647,641,703,706,685,709,704,707,686,
710,664,670,646,640
77,437,463,449,441,712,718,716,714,490,488,485,481,722,721,720,
719,711,717,715,713
*ELSET,ELSET=DEFAULT
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
73
74
75
76
77
*SURFACE,NAME=PRETENSIONSECTION
77,S1
*SURFACE,NAME=CONTACT_FACES
22,S5
20,S5
1,S5
23,S5
21,S5
14,S4
3,S5
25,S5
17,S5
33,S4
19,S6
4,S5
18,S5
32,S4
13,S6
35,S5
10,S5
24,S4
11,S5
9,S5
67,S6
66,S6
70,S6
73,S6
65,S6
69,S6
72,S6
76,S6
68,S6
71,S6
75,S6
74,S6
*SURFACE,NAME=CONTACT_FACES(2)
45,S3
60,S6
47,S3
46,S3
7,S6
49,S4
42,S3
54,S3
38,S6
55,S4
40,S3
53,S3
50,S6
39,S3
61,S3
57,S3
37,S3
59,S3
56,S3
58,S3
12,S6
8,S6
26,S6
15,S6
29,S6
16,S6
28,S6
27,S6
30,S6
34,S6
31,S6
36,S6
*SURFACE,NAME=CONTACT_FACES(4)
55,S6
40,S5
39,S5
49,S6
61,S5
42,S5
59,S5
37,S5
60,S4
47,S5
58,S5
50,S4
45,S5
7,S4
56,S5
57,S5
38,S4
46,S5
53,S5
54,S5
*SURFACE,NAME=CONTACT_FACES(3)
62,S3
64,S3
63,S3
2,S3
*SURFACE,NAME=BONDED_CONTACT_FACES(2)
41,S3
6,S3
43,S3
5,S3
*SURFACE,NAME=BONDED_CONTACT_FACES
9,S3
24,S6
11,S3
10,S3
32,S6
13,S4
35,S3
18,S3
33,S6
19,S4
4,S3
17,S3
14,S6
3,S3
25,S3
21,S3
1,S3
23,S3
20,S3
22,S3
*MATERIAL,NAME=MATERIAL
*ELASTIC,TYPE=ISOTROPIC
200000000000,0.3
*SOLID SECTION,ELSET=DEFAULT,MATERIAL=MATERIAL
*BOUNDARY
617,1,,0
617,2,,0
617,3,,0
623,1,,0
623,2,,0
623,3,,0
629,1,,0
629,2,,0
629,3,,0
635,1,,0
635,2,,0
635,3,,0
641,1,,0
641,2,,0
641,3,,0
647,1,,0
647,2,,0
647,3,,0
653,1,,0
653,2,,0
653,3,,0
659,1,,0
659,2,,0
659,3,,0
665,1,,0
665,2,,0
665,3,,0
671,1,,0
671,2,,0
671,3,,0
674,1,,0
674,2,,0
674,3,,0
677,1,,0
677,2,,0
677,3,,0
680,1,,0
680,2,,0
680,3,,0
683,1,,0
683,2,,0
683,3,,0
686,1,,0
686,2,,0
686,3,,0
689,1,,0
689,2,,0
689,3,,0
692,1,,0
692,2,,0
692,3,,0
695,1,,0
695,2,,0
695,3,,0
698,1,,0
698,2,,0
698,3,,0
701,1,,0
701,2,,0
701,3,,0
704,1,,0
704,2,,0
704,3,,0
707,1,,0
707,2,,0
707,3,,0
710,1,,0
710,2,,0
710,3,,0
*PRE-TENSION SECTION,SURFACE=PRETENSIONSECTION,NODE=723
0.000000000000E+000,1.000000000000E+000,0.000000000000E+000
*AMPLITUDE,NAME=Ax_110_1
0,0
1,-416.6666666667
*AMPLITUDE,NAME=Ax_137_2
0,0
1,-416.6666666667
*AMPLITUDE,NAME=Ax_146_3
0,0
1,-416.6666666667
*AMPLITUDE,NAME=Ax_153_4
0,0
1,-416.6666666667
*AMPLITUDE,NAME=Ax_303_5
0,0
1,1666.666666667
*AMPLITUDE,NAME=Ax_306_6
0,0
1,833.3333333333
*AMPLITUDE,NAME=Ax_309_7
0,0
1,1666.666666667
*AMPLITUDE,NAME=Ax_312_8
0,0
1,833.3333333333
*AMPLITUDE,NAME=Ax_281_9
0,0
1,-416.6666666667
*AMPLITUDE,NAME=Ax_288_A
0,0
1,-416.6666666667
*AMPLITUDE,NAME=Ax_327_B
0,0
1,833.3333333333
*AMPLITUDE,NAME=Ax_330_C
0,0
1,1666.666666667
*AMPLITUDE,NAME=Ax_333_D
0,0
1,833.3333333333
*AMPLITUDE,NAME=Ax_294_E
0,0
1,-208.3333333333
*AMPLITUDE,NAME=Ax_300_F
0,0
1,-208.3333333333
*AMPLITUDE,NAME=Ax_336_10
0,0
1,833.3333333333
*AMPLITUDE,NAME=Ax_339_11
0,0
1,833.3333333333
*AMPLITUDE,NAME=Ax_342_12
0,0
1,833.3333333333
*AMPLITUDE,NAME=Ax_159_13
0,0
1,-208.3333333333
*AMPLITUDE,NAME=Ax_197_14
0,0
1,-208.3333333333
*AMPLITUDE,NAME=Ax_316_15
0,0
1,833.3333333333
*AMPLITUDE,NAME=Ax_319_16
0,0
1,833.3333333333
*AMPLITUDE,NAME=Ax_324_17
0,0
1,833.3333333333
*CONTACT PAIR,INTERACTION=SI_18,TYPE=SURFACE TO SURFACE
CONTACT_FACES(2),CONTACT_FACES
*SURFACE INTERACTION,NAME=SI_18
*SURFACE BEHAVIOR,PRESSURE-OVERCLOSURE=LINEAR
2E+14,1
*FRICTION
0.5,2E+14
*CONTACT PAIR,INTERACTION=SI_19,TYPE=SURFACE TO SURFACE
CONTACT_FACES(3),CONTACT_FACES(4)
*SURFACE INTERACTION,NAME=SI_19
*SURFACE BEHAVIOR,PRESSURE-OVERCLOSURE=LINEAR
2E+14,1
*TIE,NAME=T_1A,POSITION TOLERANCE=0
BONDED_CONTACT_FACES(2),BONDED_CONTACT_FACES
*STEP,NLGEOM=YES,INC=100,AMPLITUDE=STEP
*STATIC,DIRECT
0.1,1,0,0
*CLOAD
723,1,30000
*CLOAD,AMPLITUDE=Ax_110_1
110,1,1
*CLOAD,AMPLITUDE=Ax_137_2
137,1,1
*CLOAD,AMPLITUDE=Ax_146_3
146,1,1
*CLOAD,AMPLITUDE=Ax_153_4
153,1,1
*CLOAD,AMPLITUDE=Ax_303_5
303,1,1
*CLOAD,AMPLITUDE=Ax_306_6
306,1,1
*CLOAD,AMPLITUDE=Ax_309_7
309,1,1
*CLOAD,AMPLITUDE=Ax_312_8
312,1,1
*CLOAD,AMPLITUDE=Ax_281_9
281,1,1
*CLOAD,AMPLITUDE=Ax_288_A
288,1,1
*CLOAD,AMPLITUDE=Ax_327_B
327,1,1
*CLOAD,AMPLITUDE=Ax_330_C
330,1,1
*CLOAD,AMPLITUDE=Ax_333_D
333,1,1
*CLOAD,AMPLITUDE=Ax_294_E
294,1,1
*CLOAD,AMPLITUDE=Ax_300_F
300,1,1
*CLOAD,AMPLITUDE=Ax_336_10
336,1,1
*CLOAD,AMPLITUDE=Ax_339_11
339,1,1
*CLOAD,AMPLITUDE=Ax_342_12
342,1,1
*CLOAD,AMPLITUDE=Ax_159_13
159,1,1
*CLOAD,AMPLITUDE=Ax_197_14
197,1,1
*CLOAD,AMPLITUDE=Ax_316_15
316,1,1
*CLOAD,AMPLITUDE=Ax_319_16
319,1,1
*CLOAD,AMPLITUDE=Ax_324_17
324,1,1
*NODE FILE,GLOBAL=YES
U,RF
*EL FILE
S,NOE
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
  {
    id: "mass",
    name: "Massepunkt",
    blurb: "*MASS + SPRINGA, f≈1.59 Hz",
    inp: `*HEADING
MASS + SPRINGA, f = sqrt(k/m)/(2π) ≈ 1.5915 Hz
*NODE
1, 0, 0, 0
2, 1, 0, 0
*ELEMENT, TYPE=SPRINGA, ELSET=S
1, 1, 2
*ELEMENT, TYPE=MASS, ELSET=M
2, 2
*SPRING, ELSET=S
100
*MASS, ELSET=M
1.0
*BOUNDARY
1, 1, 3
2, 2, 3
*STEP
*FREQUENCY
1
*NODE FILE
U
*END STEP
`,
  },
  {
    id: "gapuni",
    name: "Spalt GAPUNI",
    blurb: "geschlossener Spalt als Axialfeder, u=0.01",
    inp: `*HEADING
GAPUNI closed gap, k=1000, F=10 → u=0.01
*NODE
1, 0, 0, 0
2, 1, 0, 0
*ELEMENT, TYPE=GAPUNI, ELSET=G
1, 1, 2
*GAP, ELSET=G
0.0, 1000
*BOUNDARY
1, 1, 3
2, 2, 3
*STEP
*STATIC
*CLOAD
2, 1, 10
*NODE FILE
U, RF
*END STEP
`,
  },
  {
    id: "dashpot",
    name: "Dämpfer DASHPOTA",
    blurb: "überdämpftes SDOF, |u(T)|≪u0",
    inp: `*HEADING
DASHPOTA overdamped SDOF, ζ=1
*NODE
1, 0, 0, 0
2, 1, 0, 0
*ELEMENT, TYPE=SPRINGA, ELSET=S
1, 1, 2
*ELEMENT, TYPE=MASS, ELSET=M
2, 2
*ELEMENT, TYPE=DASHPOTA, ELSET=D
3, 1, 2
*SPRING, ELSET=S
100
*MASS, ELSET=M
1.0
*DASHPOT, ELSET=D
20
*BOUNDARY
1, 1, 3
2, 2, 3
*INITIAL CONDITIONS, TYPE=DISPLACEMENT
2, 1, 0.01
*STEP
*DYNAMIC
0.005, 0.6283185307179586
*NODE FILE
U
*END STEP
`,
  },
  {
    id: "riks",
    name: "Riks Snap-Through",
    blurb: "von Mises-Fachwerk, *STATIC, RIKS",
    inp: `*HEADING
von Mises truss — *STATIC, RIKS snap-through
*NODE
1, 0, 0, 0
2, 20, 0, 0
3, 10, 1, 0
*ELEMENT, TYPE=T3D2, ELSET=T
1, 1, 3
2, 2, 3
*MATERIAL, NAME=STEEL
*ELASTIC
210000, 0.3
*SOLID SECTION, ELSET=T, MATERIAL=STEEL
1.0
*BOUNDARY
1, 1, 3
2, 1, 3
3, 3, 3
*STEP, NLGEOM
*STATIC, RIKS
0.05, 1.0, 1e-4, 0.2, 80
*CONTROLS, MAXITER=25
*CLOAD
3, 2, -200
*NODE FILE
U, RF
*EL FILE
S
*END STEP
`,
  },
];

export const DEFAULT_EXAMPLE = EXAMPLES[0];
