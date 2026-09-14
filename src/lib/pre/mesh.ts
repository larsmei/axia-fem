import { circleLoop, pointInShape, shapeBounds, shapeHoles, shapeOutline } from "./geometry";
import type { Dim, Mesh, MeshElement, MeshNode, MeshQuality, Shape, Vec2 } from "./types";

type P = Vec2;

function subdivideLoop(loop: P[], size: number): P[] {
  if (loop.length < 2) return loop.slice();
  const out: P[] = [];
  const n = loop.length;
  for (let i = 0; i < n; i++) {
    const a = loop[i];
    const b = loop[(i + 1) % n];
    const len = Math.hypot(b.x - a.x, b.y - a.y);
    const steps = Math.max(1, Math.ceil(len / Math.max(size, 1e-6)));
    for (let s = 0; s < steps; s++) {
      const t = s / steps;
      out.push({ x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t });
    }
  }
  return out;
}

function mergeClose(pts: P[], eps: number): P[] {
  const out: P[] = [];
  for (const p of pts) {
    if (!out.some((q) => Math.hypot(q.x - p.x, q.y - p.y) < eps)) out.push(p);
  }
  return out;
}

function orient(a: P, b: P, c: P): number {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

function inCircum(a: P, b: P, c: P, p: P): boolean {
  const adx = a.x - p.x,
    ady = a.y - p.y;
  const bdx = b.x - p.x,
    bdy = b.y - p.y;
  const cdx = c.x - p.x,
    cdy = c.y - p.y;
  const det =
    (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy) -
    (bdx * bdx + bdy * bdy) * (adx * cdy - cdx * ady) +
    (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady);
  return orient(a, b, c) > 0 ? det > 1e-18 : det < -1e-18;
}

function delaunay(points: P[]): [number, number, number][] {
  const n = points.length;
  if (n < 3) return [];
  let minx = Infinity,
    miny = Infinity,
    maxx = -Infinity,
    maxy = -Infinity;
  for (const p of points) {
    minx = Math.min(minx, p.x);
    miny = Math.min(miny, p.y);
    maxx = Math.max(maxx, p.x);
    maxy = Math.max(maxy, p.y);
  }
  const dx = maxx - minx || 1;
  const dy = maxy - miny || 1;
  const d = Math.max(dx, dy) * 20;
  const mid: P = { x: (minx + maxx) / 2, y: (miny + maxy) / 2 };
  const pts = points.concat([
    { x: mid.x, y: mid.y + d },
    { x: mid.x - d, y: mid.y - d },
    { x: mid.x + d, y: mid.y - d },
  ]);
  let tris: [number, number, number][] = [[n, n + 1, n + 2]];
  for (let i = 0; i < n; i++) {
    const p = pts[i];
    const bad: [number, number, number][] = [];
    for (const t of tris) {
      if (inCircum(pts[t[0]], pts[t[1]], pts[t[2]], p)) bad.push(t);
    }
    const edges: [number, number][] = [];
    const pushE = (u: number, v: number) => {
      const ix = edges.findIndex((e) => e[0] === v && e[1] === u);
      if (ix >= 0) edges.splice(ix, 1);
      else edges.push([u, v]);
    };
    for (const t of bad) {
      pushE(t[0], t[1]);
      pushE(t[1], t[2]);
      pushE(t[2], t[0]);
    }
    const badSet = new Set(bad);
    tris = tris.filter((t) => !badSet.has(t));
    for (const [u, v] of edges) tris.push([u, v, i]);
  }
  return tris.filter((t) => t[0] < n && t[1] < n && t[2] < n);
}

function triangleCentroid(a: P, b: P, c: P): P {
  return { x: (a.x + b.x + c.x) / 3, y: (a.y + b.y + c.y) / 3 };
}

function meshPolygonShape(shape: Shape, size: number, id0: { n: number; e: number }): Mesh {
  const outline = shapeOutline(shape);
  const holes = shapeHoles(shape);
  const sz = Math.max(size, 1e-3);
  let pts = subdivideLoop(outline, sz);
  for (const h of holes) {
    const n = Math.max(12, Math.ceil((2 * Math.PI * h.r) / sz));
    pts = pts.concat(subdivideLoop(circleLoop(h.cx, h.cy, h.r, n), sz));
  }
  const b = shapeBounds(shape);
  for (let x = b.x + sz * 0.5; x < b.x + b.w; x += sz) {
    for (let y = b.y + sz * 0.5; y < b.y + b.h; y += sz) {
      const p = { x, y };
      if (pointInShape(p, shape)) pts.push(p);
    }
  }
  pts = mergeClose(pts, sz * 0.15);
  const tris = delaunay(pts);
  const nodes: MeshNode[] = pts.map((p, i) => ({ id: id0.n + i, x: p.x, y: p.y, z: 0 }));
  const elements: MeshElement[] = [];
  for (const t of tris) {
    const a = pts[t[0]],
      b = pts[t[1]],
      c = pts[t[2]];
    const cen = triangleCentroid(a, b, c);
    if (!pointInShape(cen, shape)) continue;
    if (orient(a, b, c) < 0) t.reverse();
    elements.push({
      id: ++id0.e,
      type: "CPS3",
      nodes: [id0.n + t[0], id0.n + t[1], id0.n + t[2]],
      materialId: shape.materialId,
      shapeId: shape.id,
    });
  }
  id0.n += pts.length;
  return { nodes, elements };
}

function meshRectStructured(shape: Extract<Shape, { kind: "rect" }>, size: number, id0: { n: number; e: number }): Mesh {
  const nx = Math.max(1, Math.round(shape.w / Math.max(size, 1e-6)));
  const ny = Math.max(1, Math.round(shape.h / Math.max(size, 1e-6)));
  const nodes: MeshNode[] = [];
  for (let j = 0; j <= ny; j++) {
    for (let i = 0; i <= nx; i++) {
      nodes.push({
        id: id0.n + j * (nx + 1) + i,
        x: shape.x + (shape.w * i) / nx,
        y: shape.y + (shape.h * j) / ny,
        z: 0,
      });
    }
  }
  const elements: MeshElement[] = [];
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const n00 = id0.n + j * (nx + 1) + i;
      const n10 = n00 + 1;
      const n01 = id0.n + (j + 1) * (nx + 1) + i;
      const n11 = n01 + 1;
      elements.push({
        id: ++id0.e,
        type: "CPS4",
        nodes: [n00, n10, n11, n01],
        materialId: shape.materialId,
        shapeId: shape.id,
      });
    }
  }
  id0.n += (nx + 1) * (ny + 1);
  return { nodes, elements };
}

function meshCirclePolar(shape: Extract<Shape, { kind: "circle" }>, size: number, id0: { n: number; e: number }): Mesh {
  const r = Math.max(shape.r, 1e-6);
  const nr = Math.max(2, Math.round(r / Math.max(size, 1e-6)));
  const nt = Math.max(8, Math.round((2 * Math.PI * r) / Math.max(size, 1e-6)));
  const nodes: MeshNode[] = [{ id: id0.n, x: shape.cx, y: shape.cy, z: 0 }];
  for (let k = 1; k <= nr; k++) {
    const rr = (r * k) / nr;
    for (let t = 0; t < nt; t++) {
      const a = (t / nt) * Math.PI * 2;
      nodes.push({
        id: id0.n + 1 + (k - 1) * nt + t,
        x: shape.cx + rr * Math.cos(a),
        y: shape.cy + rr * Math.sin(a),
        z: 0,
      });
    }
  }
  const elements: MeshElement[] = [];
  const ring = (k: number, t: number) => id0.n + 1 + (k - 1) * nt + (t % nt);
  for (let t = 0; t < nt; t++) {
    elements.push({
      id: ++id0.e,
      type: "CPS3",
      nodes: [id0.n, ring(1, t), ring(1, t + 1)],
      materialId: shape.materialId,
      shapeId: shape.id,
    });
  }
  for (let k = 1; k < nr; k++) {
    for (let t = 0; t < nt; t++) {
      elements.push({
        id: ++id0.e,
        type: "CPS4",
        nodes: [ring(k, t), ring(k + 1, t), ring(k + 1, t + 1), ring(k, t + 1)],
        materialId: shape.materialId,
        shapeId: shape.id,
      });
    }
  }
  id0.n += 1 + nr * nt;
  return { nodes, elements };
}

function meshBoxStructured(
  x: number,
  y: number,
  z: number,
  w: number,
  h: number,
  d: number,
  size: number,
  materialId: string,
  shapeId: string,
  id0: { n: number; e: number },
): Mesh {
  const nx = Math.max(1, Math.round(w / Math.max(size, 1e-6)));
  const ny = Math.max(1, Math.round(h / Math.max(size, 1e-6)));
  const nz = Math.max(1, Math.round(d / Math.max(size, 1e-6)));
  const sx = nx + 1,
    sy = ny + 1;
  const nid = (i: number, j: number, k: number) => id0.n + k * sx * sy + j * sx + i;
  const nodes: MeshNode[] = [];
  for (let k = 0; k <= nz; k++) {
    for (let j = 0; j <= ny; j++) {
      for (let i = 0; i <= nx; i++) {
        nodes.push({
          id: nid(i, j, k),
          x: x + (w * i) / nx,
          y: y + (h * j) / ny,
          z: z + (d * k) / nz,
        });
      }
    }
  }
  const elements: MeshElement[] = [];
  for (let k = 0; k < nz; k++) {
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        elements.push({
          id: ++id0.e,
          type: "C3D8",
          nodes: [
            nid(i, j, k),
            nid(i + 1, j, k),
            nid(i + 1, j + 1, k),
            nid(i, j + 1, k),
            nid(i, j, k + 1),
            nid(i + 1, j, k + 1),
            nid(i + 1, j + 1, k + 1),
            nid(i, j + 1, k + 1),
          ],
          materialId,
          shapeId,
        });
      }
    }
  }
  id0.n += sx * sy * (nz + 1);
  return { nodes, elements };
}

function meshCylinderStructured(
  cx: number,
  cy: number,
  cz: number,
  r: number,
  height: number,
  axis: "x" | "y" | "z",
  size: number,
  materialId: string,
  shapeId: string,
  id0: { n: number; e: number },
): Mesh {
  const rr = Math.max(r, 1e-6);
  const nr = Math.max(2, Math.round(rr / Math.max(size, 1e-6)));
  const nt = Math.max(8, Math.round((2 * Math.PI * rr) / Math.max(size, 1e-6)));
  const nz = Math.max(1, Math.round(height / Math.max(size, 1e-6)));
  const hh = height / 2;
  const map = (px: number, py: number, pz: number) => {
    if (axis === "z") return { x: cx + px, y: cy + py, z: cz + pz };
    if (axis === "y") return { x: cx + px, y: cy + pz, z: cz + py };
    return { x: cx + pz, y: cy + px, z: cz + py };
  };
  const centerOf = (k: number) => id0.n + k;
  const ringOf = (j: number, t: number, k: number) => id0.n + (nz + 1) + k * nr * nt + (j - 1) * nt + (t % nt);
  const nodes: MeshNode[] = [];
  for (let k = 0; k <= nz; k++) {
    const z = -hh + (height * k) / nz;
    const p = map(0, 0, z);
    nodes.push({ id: centerOf(k), x: p.x, y: p.y, z: p.z });
  }
  for (let k = 0; k <= nz; k++) {
    const z = -hh + (height * k) / nz;
    for (let j = 1; j <= nr; j++) {
      const rad = (rr * j) / nr;
      for (let t = 0; t < nt; t++) {
        const a = (t / nt) * Math.PI * 2;
        const p = map(rad * Math.cos(a), rad * Math.sin(a), z);
        nodes.push({ id: ringOf(j, t, k), x: p.x, y: p.y, z: p.z });
      }
    }
  }
  const elements: MeshElement[] = [];
  for (let k = 0; k < nz; k++) {
    for (let t = 0; t < nt; t++) {
      elements.push({
        id: ++id0.e,
        type: "C3D6",
        nodes: [centerOf(k), ringOf(1, t, k), ringOf(1, t + 1, k), centerOf(k + 1), ringOf(1, t, k + 1), ringOf(1, t + 1, k + 1)],
        materialId,
        shapeId,
      });
    }
    for (let j = 1; j < nr; j++) {
      for (let t = 0; t < nt; t++) {
        elements.push({
          id: ++id0.e,
          type: "C3D8",
          nodes: [
            ringOf(j, t, k),
            ringOf(j + 1, t, k),
            ringOf(j + 1, t + 1, k),
            ringOf(j, t + 1, k),
            ringOf(j, t, k + 1),
            ringOf(j + 1, t, k + 1),
            ringOf(j + 1, t + 1, k + 1),
            ringOf(j, t + 1, k + 1),
          ],
          materialId,
          shapeId,
        });
      }
    }
  }
  id0.n += nz + 1 + (nz + 1) * nr * nt;
  return { nodes, elements };
}

function cubeToBall(u: number, v: number, w: number, r: number): { x: number; y: number; z: number } {
  const inf = Math.max(Math.abs(u), Math.abs(v), Math.abs(w));
  const l2 = Math.hypot(u, v, w);
  if (l2 < 1e-18) return { x: 0, y: 0, z: 0 };
  const s = (r * inf) / l2;
  return { x: u * s, y: v * s, z: w * s };
}

function meshSphereHex(
  cx: number,
  cy: number,
  cz: number,
  r: number,
  size: number,
  materialId: string,
  shapeId: string,
  id0: { n: number; e: number },
): Mesh {
  const n = Math.max(2, Math.round((2 * r) / Math.max(size, 1e-6)));
  const nid = (i: number, j: number, k: number) => id0.n + k * (n + 1) * (n + 1) + j * (n + 1) + i;
  const nodes: MeshNode[] = [];
  for (let k = 0; k <= n; k++) {
    for (let j = 0; j <= n; j++) {
      for (let i = 0; i <= n; i++) {
        const u = -1 + (2 * i) / n;
        const v = -1 + (2 * j) / n;
        const w = -1 + (2 * k) / n;
        const p = cubeToBall(u, v, w, r);
        nodes.push({ id: nid(i, j, k), x: cx + p.x, y: cy + p.y, z: cz + p.z });
      }
    }
  }
  const elements: MeshElement[] = [];
  for (let k = 0; k < n; k++) {
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < n; i++) {
        elements.push({
          id: ++id0.e,
          type: "C3D8",
          nodes: [
            nid(i, j, k),
            nid(i + 1, j, k),
            nid(i + 1, j + 1, k),
            nid(i, j + 1, k),
            nid(i, j, k + 1),
            nid(i + 1, j, k + 1),
            nid(i + 1, j + 1, k + 1),
            nid(i, j + 1, k + 1),
          ],
          materialId,
          shapeId,
        });
      }
    }
  }
  id0.n += (n + 1) ** 3;
  return { nodes, elements };
}

function extrude2d(base: Mesh, depth: number, size: number, id0: { n: number; e: number }): Mesh {
  const nz = Math.max(1, Math.round(depth / Math.max(size, 1e-6)));
  const dz = depth / nz;
  const nBase = base.nodes.length;
  const indexOf = new Map(base.nodes.map((n, i) => [n.id, i]));
  const nodes: MeshNode[] = [];
  for (let k = 0; k <= nz; k++) {
    const z = k * dz;
    for (const n of base.nodes) {
      nodes.push({ id: id0.n + k * nBase + indexOf.get(n.id)!, x: n.x, y: n.y, z });
    }
  }
  const nid = (orig: number, k: number) => id0.n + k * nBase + indexOf.get(orig)!;
  const elements: MeshElement[] = [];
  for (let k = 0; k < nz; k++) {
    for (const el of base.elements) {
      if (el.type === "CPS4" && el.nodes.length >= 4) {
        const [a, b, c, d] = el.nodes;
        elements.push({
          id: ++id0.e,
          type: "C3D8",
          nodes: [nid(a, k), nid(b, k), nid(c, k), nid(d, k), nid(a, k + 1), nid(b, k + 1), nid(c, k + 1), nid(d, k + 1)],
          materialId: el.materialId,
          shapeId: el.shapeId,
        });
      } else if (el.nodes.length >= 3) {
        const [a, b, c] = el.nodes;
        elements.push({
          id: ++id0.e,
          type: "C3D6",
          nodes: [nid(a, k), nid(b, k), nid(c, k), nid(a, k + 1), nid(b, k + 1), nid(c, k + 1)],
          materialId: el.materialId,
          shapeId: el.shapeId,
        });
      }
    }
  }
  id0.n += nBase * (nz + 1);
  return { nodes, elements };
}

function meshOne(shape: Shape, size: number, dim: Dim, id0: { n: number; e: number }): Mesh {
  if (shape.kind === "box") {
    return meshBoxStructured(shape.x, shape.y, shape.z, shape.w, shape.h, shape.d, size, shape.materialId, shape.id, id0);
  }
  if (shape.kind === "cylinder") {
    return meshCylinderStructured(
      shape.cx,
      shape.cy,
      shape.cz,
      shape.r,
      shape.height,
      shape.axis,
      size,
      shape.materialId,
      shape.id,
      id0,
    );
  }
  if (shape.kind === "sphere") {
    return meshSphereHex(shape.cx, shape.cy, shape.cz, shape.r, size, shape.materialId, shape.id, id0);
  }
  if (dim === "3d") {
    const depth = Math.max(shape.depth ?? 1, size * 0.5);
    if (shape.kind === "rect" && (!shape.holes || shape.holes.length === 0)) {
      return meshBoxStructured(shape.x, shape.y, 0, shape.w, shape.h, depth, size, shape.materialId, shape.id, id0);
    }
    if (shape.kind === "circle") {
      return meshCylinderStructured(shape.cx, shape.cy, depth / 2, shape.r, depth, "z", size, shape.materialId, shape.id, id0);
    }
    const flat = meshPolygonShape(shape, size, { n: 1, e: 0 });
    return extrude2d(flat, depth, size, id0);
  }
  if (shape.kind === "rect" && (!shape.holes || shape.holes.length === 0)) return meshRectStructured(shape, size, id0);
  if (shape.kind === "circle") return meshCirclePolar(shape, size, id0);
  return meshPolygonShape(shape, size, id0);
}

export function meshShapes(shapes: Shape[], size: number, dim: Dim = "2d"): Mesh {
  const id0 = { n: 1, e: 0 };
  const nodes: MeshNode[] = [];
  const elements: MeshElement[] = [];
  for (const s of shapes) {
    const part = meshOne(s, size, dim, id0);
    nodes.push(...part.nodes);
    elements.push(...part.elements);
  }
  return { nodes, elements };
}

function cornerAngles(pts: { x: number; y: number; z?: number }[]): number[] {
  const n = pts.length;
  const out: number[] = [];
  for (let i = 0; i < n; i++) {
    const a = pts[(i + n - 1) % n];
    const b = pts[i];
    const c = pts[(i + 1) % n];
    const v1x = a.x - b.x,
      v1y = a.y - b.y,
      v1z = (a.z ?? 0) - (b.z ?? 0);
    const v2x = c.x - b.x,
      v2y = c.y - b.y,
      v2z = (c.z ?? 0) - (b.z ?? 0);
    const d1 = Math.hypot(v1x, v1y, v1z);
    const d2 = Math.hypot(v2x, v2y, v2z);
    if (d1 < 1e-12 || d2 < 1e-12) continue;
    const cos = Math.max(-1, Math.min(1, (v1x * v2x + v1y * v2y + v1z * v2z) / (d1 * d2)));
    out.push((Math.acos(cos) * 180) / Math.PI);
  }
  return out;
}

const HEX_FACES = [
  [0, 1, 2, 3],
  [4, 5, 6, 7],
  [0, 1, 5, 4],
  [1, 2, 6, 5],
  [2, 3, 7, 6],
  [3, 0, 4, 7],
];

export function meshQuality(mesh: Mesh): MeshQuality {
  const byId = new Map(mesh.nodes.map((n) => [n.id, n]));
  let minAngle = 180;
  let maxAspect = 1;
  let nBad = 0;
  const types = [...new Set(mesh.elements.map((e) => e.type))];
  for (const el of mesh.elements) {
    const pts = el.nodes.map((id) => byId.get(id)).filter((n): n is MeshNode => !!n);
    if (pts.length < 3) continue;
    const faces: MeshNode[][] = [];
    if (el.type === "C3D8" && pts.length >= 8) {
      for (const f of HEX_FACES) faces.push(f.map((i) => pts[i]));
    } else if (el.type === "C3D6" && pts.length >= 6) {
      faces.push([pts[0], pts[1], pts[2]], [pts[3], pts[4], pts[5]]);
    } else {
      faces.push(pts);
    }
    const edges: number[] = [];
    for (let i = 0; i < pts.length; i++) {
      for (let j = i + 1; j < pts.length; j++) {
        const a = pts[i],
          b = pts[j];
        const d = Math.hypot(b.x - a.x, b.y - a.y, (b.z ?? 0) - (a.z ?? 0));
        if (d > 1e-12) edges.push(d);
      }
    }
    const aspect = edges.length ? Math.max(...edges) / Math.max(Math.min(...edges), 1e-12) : 1;
    let amin = 180;
    for (const f of faces) {
      const angs = cornerAngles(f);
      if (angs.length) amin = Math.min(amin, ...angs);
    }
    minAngle = Math.min(minAngle, amin);
    maxAspect = Math.max(maxAspect, aspect);
    if (amin < 12 || aspect > 12) nBad += 1;
  }
  return {
    nnode: mesh.nodes.length,
    nelem: mesh.elements.length,
    minAngle,
    maxAspect,
    nBad,
    types,
  };
}

export function elemTypeSummary(mesh: Mesh | null): string {
  if (!mesh?.elements.length) return "—";
  const q = meshQuality(mesh);
  return q.types.join(" / ");
}
