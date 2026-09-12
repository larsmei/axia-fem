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

export const EXAMPLES: Example[] = [
  {
    id: "patch",
    name: "Zugstab C3D8",
    blurb: "Patch-Test, ein Hexaeder",
    inp: tensionPatch.trimStart(),
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
];

export const DEFAULT_EXAMPLE = EXAMPLES[2];
