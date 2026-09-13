import type { EdgeName, Hole, Mesh, MeshNode, Shape, Vec2 } from "./types";

export function dist(a: Vec2, b: Vec2): number {
  const dx = a.x - b.x;
  const dy = a.y - b.y;
  return Math.hypot(dx, dy);
}

export function lerp(a: Vec2, b: Vec2, t: number): Vec2 {
  return { x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t };
}

export type Bounds = { x: number; y: number; w: number; h: number };

export function shapeBounds(s: Shape): Bounds {
  if (s.kind === "rect") return { x: s.x, y: s.y, w: s.w, h: s.h };
  if (s.kind === "circle") {
    return { x: s.cx - s.r, y: s.cy - s.r, w: s.r * 2, h: s.r * 2 };
  }
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
  if (!s.points.length) return { x: 0, y: 0, w: 0, h: 0 };
  return { x: minx, y: miny, w: maxx - minx, h: maxy - miny };
}

export function modelBounds(shapes: Shape[]): Bounds {
  if (!shapes.length) return { x: -20, y: -20, w: 140, h: 80 };
  let minx = Infinity,
    miny = Infinity,
    maxx = -Infinity,
    maxy = -Infinity;
  for (const s of shapes) {
    const b = shapeBounds(s);
    minx = Math.min(minx, b.x);
    miny = Math.min(miny, b.y);
    maxx = Math.max(maxx, b.x + b.w);
    maxy = Math.max(maxy, b.y + b.h);
  }
  const pad = Math.max(10, 0.08 * Math.max(maxx - minx, maxy - miny));
  return { x: minx - pad, y: miny - pad, w: maxx - minx + 2 * pad, h: maxy - miny + 2 * pad };
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
  if (s.kind === "rect") {
    return [
      { x: s.x, y: s.y },
      { x: s.x + s.w, y: s.y },
      { x: s.x + s.w, y: s.y + s.h },
      { x: s.x, y: s.y + s.h },
    ];
  }
  if (s.kind === "circle") {
    const n = Math.max(24, Math.ceil((2 * Math.PI * s.r) / Math.max(s.r / 8, 1)));
    return circleLoop(s.cx, s.cy, s.r, n);
  }
  return s.points;
}

export function shapeHoles(s: Shape): Hole[] {
  if (s.kind === "circle") return [];
  return s.holes ?? [];
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
  if (s.kind === "rect") {
    if (p.x < s.x || p.y < s.y || p.x > s.x + s.w || p.y > s.y + s.h) return false;
    for (const h of s.holes) if (dist(p, { x: h.cx, y: h.cy }) <= h.r) return false;
    return true;
  }
  if (s.kind === "circle") return dist(p, { x: s.cx, y: s.cy }) <= s.r;
  if (!pointInPolygon(p, s.points)) return false;
  for (const h of s.holes) if (dist(p, { x: h.cx, y: h.cy }) <= h.r) return false;
  return true;
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

export function edgeNodes(mesh: Mesh, shape: Shape, edge: EdgeName, tol: number): MeshNode[] {
  const b = shapeBounds(shape);
  const outline = shapeOutline(shape);
  return mesh.nodes.filter((n) => {
    if (!nodeOnShape(n, shape, tol * 2)) return false;
    if (edge === "boundary") return true;
    if (edge === "left") return Math.abs(n.x - b.x) <= tol;
    if (edge === "right") return Math.abs(n.x - (b.x + b.w)) <= tol;
    if (edge === "bottom") return Math.abs(n.y - b.y) <= tol;
    if (edge === "top") return Math.abs(n.y - (b.y + b.h)) <= tol;
    return outline.some((_, i) => distToSegment(n, outline[i], outline[(i + 1) % outline.length]) <= tol);
  });
}

function nodeOnShape(n: MeshNode, shape: Shape, tol: number): boolean {
  if (pointInShape(n, shape)) return true;
  const outline = shapeOutline(shape);
  for (let i = 0; i < outline.length; i++) {
    if (distToSegment(n, outline[i], outline[(i + 1) % outline.length]) <= tol) return true;
  }
  return false;
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
