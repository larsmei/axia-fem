import type { EdgeName, FaceName, Hole, Mesh, MeshNode, Shape, Vec2, Vec3 } from "./types";

export function dist(a: Vec2, b: Vec2): number {
  return Math.hypot(a.x - b.x, a.y - b.y);
}

export function dist3(a: Vec3, b: Vec3): number {
  return Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z);
}

export function lerp(a: Vec2, b: Vec2, t: number): Vec2 {
  return { x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t };
}

export type Bounds = { x: number; y: number; w: number; h: number };
export type Bounds3 = { x: number; y: number; z: number; w: number; h: number; d: number };

export function isSolid3(s: Shape): boolean {
  return s.kind === "box" || s.kind === "cylinder" || s.kind === "sphere";
}

export function shapeDepth(s: Shape): number {
  if (s.kind === "box") return s.d;
  if (s.kind === "cylinder") return s.height;
  if (s.kind === "sphere") return s.r * 2;
  return s.depth ?? 0;
}

export function shapeBounds(s: Shape): Bounds {
  const b = shapeBounds3(s);
  return { x: b.x, y: b.y, w: b.w, h: b.h };
}

export function shapeBounds3(s: Shape): Bounds3 {
  if (s.kind === "rect") return { x: s.x, y: s.y, z: 0, w: s.w, h: s.h, d: s.depth ?? 0 };
  if (s.kind === "circle") {
    return { x: s.cx - s.r, y: s.cy - s.r, z: 0, w: s.r * 2, h: s.r * 2, d: s.depth ?? 0 };
  }
  if (s.kind === "polygon") {
    let minx = Infinity,
      miny = Infinity,
      maxx = -Infinity,
      maxy = -Infinity;
    for (const p of s.points) {
      minx = Math.min(minx, p.x);
      miny = Math.min(miny, p.y);
      maxx = Math.max(maxx, p.x);
      maxy = Math.max(maxy, p.y);
    }
    if (!s.points.length) return { x: 0, y: 0, z: 0, w: 0, h: 0, d: s.depth ?? 0 };
    return { x: minx, y: miny, z: 0, w: maxx - minx, h: maxy - miny, d: s.depth ?? 0 };
  }
  if (s.kind === "box") return { x: s.x, y: s.y, z: s.z, w: s.w, h: s.h, d: s.d };
  if (s.kind === "cylinder") {
    const r = s.r;
    const hh = s.height / 2;
    if (s.axis === "x") return { x: s.cx - hh, y: s.cy - r, z: s.cz - r, w: s.height, h: r * 2, d: r * 2 };
    if (s.axis === "y") return { x: s.cx - r, y: s.cy - hh, z: s.cz - r, w: r * 2, h: s.height, d: r * 2 };
    return { x: s.cx - r, y: s.cy - r, z: s.cz - hh, w: r * 2, h: r * 2, d: s.height };
  }
  return { x: s.cx - s.r, y: s.cy - s.r, z: s.cz - s.r, w: s.r * 2, h: s.r * 2, d: s.r * 2 };
}

export function modelBounds(shapes: Shape[]): Bounds {
  const b = modelBounds3(shapes);
  return { x: b.x, y: b.y, w: b.w, h: b.h };
}

export function modelBounds3(shapes: Shape[]): Bounds3 {
  if (!shapes.length) return { x: -20, y: -20, z: -10, w: 140, h: 80, d: 40 };
  let minx = Infinity,
    miny = Infinity,
    minz = Infinity,
    maxx = -Infinity,
    maxy = -Infinity,
    maxz = -Infinity;
  for (const s of shapes) {
    const b = shapeBounds3(s);
    minx = Math.min(minx, b.x);
    miny = Math.min(miny, b.y);
    minz = Math.min(minz, b.z);
    maxx = Math.max(maxx, b.x + b.w);
    maxy = Math.max(maxy, b.y + b.h);
    maxz = Math.max(maxz, b.z + b.d);
  }
  const span = Math.max(maxx - minx, maxy - miny, maxz - minz, 1);
  const pad = Math.max(10, 0.08 * span);
  return {
    x: minx - pad,
    y: miny - pad,
    z: minz - pad,
    w: maxx - minx + 2 * pad,
    h: maxy - miny + 2 * pad,
    d: maxz - minz + 2 * pad,
  };
}

export function circleLoop(cx: number, cy: number, r: number, n = 32): Vec2[] {
  const pts: Vec2[] = [];
  const k = Math.max(8, n);
  for (let i = 0; i < k; i++) {
    const a = (i / k) * Math.PI * 2;
    pts.push({ x: cx + r * Math.cos(a), y: cy + r * Math.sin(a) });
  }
  return pts;
}

export function shapeOutline(s: Shape): Vec2[] {
  if (s.kind === "rect" || s.kind === "box") {
    const x = s.x,
      y = s.y,
      w = s.kind === "rect" ? s.w : s.w,
      h = s.kind === "rect" ? s.h : s.h;
    return [
      { x, y },
      { x: x + w, y },
      { x: x + w, y: y + h },
      { x, y: y + h },
    ];
  }
  if (s.kind === "circle") {
    const n = Math.max(24, Math.ceil((2 * Math.PI * s.r) / Math.max(s.r / 8, 1)));
    return circleLoop(s.cx, s.cy, s.r, n);
  }
  if (s.kind === "cylinder") {
    const n = Math.max(24, Math.ceil((2 * Math.PI * s.r) / Math.max(s.r / 8, 1)));
    if (s.axis === "z") return circleLoop(s.cx, s.cy, s.r, n);
    if (s.axis === "y") return circleLoop(s.cx, s.cz, s.r, n);
    return circleLoop(s.cy, s.cz, s.r, n);
  }
  if (s.kind === "sphere") {
    const n = Math.max(24, Math.ceil((2 * Math.PI * s.r) / Math.max(s.r / 8, 1)));
    return circleLoop(s.cx, s.cy, s.r, n);
  }
  return s.points;
}

export function shapeHoles(s: Shape): Hole[] {
  if (s.kind === "rect" || s.kind === "polygon") return s.holes ?? [];
  return [];
}

export function pointInPolygon(p: Vec2, ring: Vec2[]): boolean {
  let inside = false;
  for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
    const a = ring[i];
    const b = ring[j];
    const hit = a.y > p.y !== b.y > p.y && p.x < ((b.x - a.x) * (p.y - a.y)) / (b.y - a.y + 1e-18) + a.x;
    if (hit) inside = !inside;
  }
  return inside;
}

export function pointInShape(p: Vec2, s: Shape): boolean {
  if (s.kind === "rect" || s.kind === "box") {
    const w = s.w;
    const h = s.h;
    if (p.x < s.x || p.y < s.y || p.x > s.x + w || p.y > s.y + h) return false;
    if (s.kind === "rect") for (const hole of s.holes) if (dist(p, { x: hole.cx, y: hole.cy }) <= hole.r) return false;
    return true;
  }
  if (s.kind === "circle") return dist(p, { x: s.cx, y: s.cy }) <= s.r;
  if (s.kind === "cylinder" && s.axis === "z") return dist(p, { x: s.cx, y: s.cy }) <= s.r;
  if (s.kind === "sphere") return dist(p, { x: s.cx, y: s.cy }) <= s.r;
  if (s.kind === "polygon") {
    if (!pointInPolygon(p, s.points)) return false;
    for (const hole of s.holes) if (dist(p, { x: hole.cx, y: hole.cy }) <= hole.r) return false;
    return true;
  }
  return false;
}

export function hitShape(shapes: Shape[], p: Vec2, slop: number): Shape | null {
  for (let i = shapes.length - 1; i >= 0; i--) {
    const s = shapes[i];
    if (pointInShape(p, s)) return s;
    const outline = shapeOutline(s);
    for (let k = 0; k < outline.length; k++) {
      const a = outline[k];
      const b = outline[(k + 1) % outline.length];
      if (distToSegment(p, a, b) <= slop) return s;
    }
  }
  return null;
}

function distToSegment(p: Vec2, a: Vec2, b: Vec2): number {
  const vx = b.x - a.x;
  const vy = b.y - a.y;
  const l2 = vx * vx + vy * vy;
  if (l2 < 1e-18) return dist(p, a);
  let t = ((p.x - a.x) * vx + (p.y - a.y) * vy) / l2;
  t = Math.max(0, Math.min(1, t));
  return dist(p, { x: a.x + t * vx, y: a.y + t * vy });
}

export function snap(v: number, step: number): number {
  if (step <= 0) return v;
  return Math.round(v / step) * step;
}

export function snapPt(p: Vec2, step: number): Vec2 {
  return { x: snap(p.x, step), y: snap(p.y, step) };
}

export function snapPt3(p: Vec3, step: number): Vec3 {
  return { x: snap(p.x, step), y: snap(p.y, step), z: snap(p.z, step) };
}

export function mapEdgeToFace(edge: EdgeName): FaceName {
  if (edge === "left") return "xmin";
  if (edge === "right") return "xmax";
  if (edge === "bottom") return "ymin";
  if (edge === "top") return "ymax";
  return "lateral";
}

export function faceLabel(face: FaceName | EdgeName): string {
  const m: Record<string, string> = {
    left: "X−",
    right: "X+",
    bottom: "Y−",
    top: "Y+",
    boundary: "Rand",
    xmin: "X−",
    xmax: "X+",
    ymin: "Y−",
    ymax: "Y+",
    zmin: "Z−",
    zmax: "Z+",
    lateral: "Mantel",
  };
  return m[face] ?? face;
}

function nodeOnShape(n: MeshNode, shape: Shape, tol: number): boolean {
  if (pointInShape(n, shape)) return true;
  const outline = shapeOutline(shape);
  for (let i = 0; i < outline.length; i++) {
    if (distToSegment(n, outline[i], outline[(i + 1) % outline.length]) <= tol) return true;
  }
  return false;
}

export function edgeNodes(mesh: Mesh, shape: Shape, edge: EdgeName, tol: number): MeshNode[] {
  return faceNodes(mesh, shape, mapEdgeToFace(edge), tol);
}

export function faceNodes(mesh: Mesh, shape: Shape, face: FaceName, tol: number): MeshNode[] {
  const b = shapeBounds3(shape);
  const t = Math.max(tol, 1e-4);
  return mesh.nodes.filter((n) => {
    if (shape.kind === "sphere") {
      const d = dist3(n, { x: shape.cx, y: shape.cy, z: shape.cz });
      if (Math.abs(d - shape.r) > t * 4 && d < shape.r - t) return false;
    } else if (shape.kind === "cylinder") {
      // accept nodes belonging to this body via bbox
      if (n.x < b.x - t || n.y < b.y - t || n.z < b.z - t || n.x > b.x + b.w + t || n.y > b.y + b.h + t || n.z > b.z + b.d + t)
        return false;
    } else if (!nodeOnShape(n, shape, t * 2) && (n.z < b.z - t || n.z > b.z + b.d + t)) {
      if (n.x < b.x - t || n.y < b.y - t || n.x > b.x + b.w + t || n.y > b.y + b.h + t) return false;
    }

    const inBbox =
      n.x >= b.x - t &&
      n.y >= b.y - t &&
      n.z >= b.z - t &&
      n.x <= b.x + b.w + t &&
      n.y <= b.y + b.h + t &&
      n.z <= b.z + b.d + t;
    if (!inBbox && shape.kind !== "sphere") return false;

    if (face === "xmin") return Math.abs(n.x - b.x) <= t;
    if (face === "xmax") return Math.abs(n.x - (b.x + b.w)) <= t;
    if (face === "ymin") return Math.abs(n.y - b.y) <= t;
    if (face === "ymax") return Math.abs(n.y - (b.y + b.h)) <= t;
    if (face === "zmin") return Math.abs(n.z - b.z) <= t;
    if (face === "zmax") return Math.abs(n.z - (b.z + b.d)) <= t;
    if (face === "lateral") {
      if (shape.kind === "cylinder") {
        const ax = shape.axis;
        const rr = shape.r;
        if (ax === "z") return Math.abs(Math.hypot(n.x - shape.cx, n.y - shape.cy) - rr) <= t * 2;
        if (ax === "y") return Math.abs(Math.hypot(n.x - shape.cx, n.z - shape.cz) - rr) <= t * 2;
        return Math.abs(Math.hypot(n.y - shape.cy, n.z - shape.cz) - rr) <= t * 2;
      }
      if (shape.kind === "circle") return Math.abs(Math.hypot(n.x - shape.cx, n.y - shape.cy) - shape.r) <= t * 2;
      if (shape.kind === "sphere") return Math.abs(dist3(n, { x: shape.cx, y: shape.cy, z: shape.cz }) - shape.r) <= t * 3;
      return true;
    }
    return false;
  });
}

export type FaceQuad = { face: FaceName; corners: Vec3[]; center: Vec3; normal: Vec3 };

export function shapeFaces(s: Shape): FaceQuad[] {
  const b = shapeBounds3(s);
  const { x, y, z, w, h, d } = b;
  const quad = (face: FaceName, c: Vec3[], n: Vec3): FaceQuad => ({
    face,
    corners: c,
    center: {
      x: (c[0].x + c[1].x + c[2].x + c[3].x) / 4,
      y: (c[0].y + c[1].y + c[2].y + c[3].y) / 4,
      z: (c[0].z + c[1].z + c[2].z + c[3].z) / 4,
    },
    normal: n,
  });
  if (s.kind === "cylinder") {
    const hh = s.height / 2;
    const r = s.r;
    if (s.axis === "z") {
      return [
        quad(
          "zmin",
          [
            { x: s.cx - r, y: s.cy - r, z: s.cz - hh },
            { x: s.cx + r, y: s.cy - r, z: s.cz - hh },
            { x: s.cx + r, y: s.cy + r, z: s.cz - hh },
            { x: s.cx - r, y: s.cy + r, z: s.cz - hh },
          ],
          { x: 0, y: 0, z: -1 },
        ),
        quad(
          "zmax",
          [
            { x: s.cx - r, y: s.cy - r, z: s.cz + hh },
            { x: s.cx + r, y: s.cy - r, z: s.cz + hh },
            { x: s.cx + r, y: s.cy + r, z: s.cz + hh },
            { x: s.cx - r, y: s.cy + r, z: s.cz + hh },
          ],
          { x: 0, y: 0, z: 1 },
        ),
      ];
    }
  }
  return [
    quad(
      "xmin",
      [
        { x, y, z },
        { x, y: y + h, z },
        { x, y: y + h, z: z + d },
        { x, y, z: z + d },
      ],
      { x: -1, y: 0, z: 0 },
    ),
    quad(
      "xmax",
      [
        { x: x + w, y, z },
        { x: x + w, y: y + h, z },
        { x: x + w, y: y + h, z: z + d },
        { x: x + w, y, z: z + d },
      ],
      { x: 1, y: 0, z: 0 },
    ),
    quad(
      "ymin",
      [
        { x, y, z },
        { x: x + w, y, z },
        { x: x + w, y, z: z + d },
        { x, y, z: z + d },
      ],
      { x: 0, y: -1, z: 0 },
    ),
    quad(
      "ymax",
      [
        { x, y: y + h, z },
        { x: x + w, y: y + h, z },
        { x: x + w, y: y + h, z: z + d },
        { x, y: y + h, z: z + d },
      ],
      { x: 0, y: 1, z: 0 },
    ),
    quad(
      "zmin",
      [
        { x, y, z },
        { x: x + w, y, z },
        { x: x + w, y: y + h, z },
        { x, y: y + h, z },
      ],
      { x: 0, y: 0, z: -1 },
    ),
    quad(
      "zmax",
      [
        { x, y, z: z + d },
        { x: x + w, y, z: z + d },
        { x: x + w, y: y + h, z: z + d },
        { x, y: y + h, z: z + d },
      ],
      { x: 0, y: 0, z: 1 },
    ),
  ];
}

export function nearestNode(mesh: Mesh, p: Vec2, maxDist: number): MeshNode | null {
  let best: MeshNode | null = null;
  let d0 = maxDist;
  for (const n of mesh.nodes) {
    const d = dist(n, p);
    if (d < d0) {
      d0 = d;
      best = n;
    }
  }
  return best;
}

export function nearestNode3(mesh: Mesh, origin: Vec3, dir: Vec3, maxDist: number): MeshNode | null {
  let best: MeshNode | null = null;
  let bestD = maxDist;
  const dl = Math.hypot(dir.x, dir.y, dir.z) || 1;
  const dx = dir.x / dl,
    dy = dir.y / dl,
    dz = dir.z / dl;
  for (const n of mesh.nodes) {
    const vx = n.x - origin.x,
      vy = n.y - origin.y,
      vz = n.z - origin.z;
    const t = vx * dx + vy * dy + vz * dz;
    if (t < 0) continue;
    const px = origin.x + t * dx,
      py = origin.y + t * dy,
      pz = origin.z + t * dz;
    const d = Math.hypot(n.x - px, n.y - py, n.z - pz);
    if (d < bestD) {
      bestD = d;
      best = n;
    }
  }
  return best;
}

export function shapeKindLabel(s: Shape): string {
  if (s.kind === "rect") return "Rechteck";
  if (s.kind === "circle") return "Kreis";
  if (s.kind === "polygon") return "Polygon";
  if (s.kind === "box") return "Quader";
  if (s.kind === "cylinder") return "Zylinder";
  return "Kugel";
}
