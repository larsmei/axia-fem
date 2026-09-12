import { useEffect, useRef } from "react";
import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { sampleColor } from "@/lib/colormap";
import type { FemResult } from "@/lib/solver";

export type FieldId = "vm" | "u" | "ux" | "uy" | "uz" | "sxx" | "syy" | "szz" | "sxy";

type Props = {
  mesh: FemResult | null;
  result: FemResult | null;
  field: FieldId;
  deformed: boolean;
  scale: number;
};

const HEX20_FACES = [
  [0, 1, 2, 3, 8, 9, 10, 11],
  [4, 7, 6, 5, 15, 14, 13, 12],
  [0, 4, 5, 1, 16, 12, 17, 8],
  [1, 5, 6, 2, 17, 13, 18, 9],
  [2, 6, 7, 3, 18, 14, 19, 10],
  [3, 7, 4, 0, 19, 15, 16, 11],
];
const TET10_FACES = [
  [0, 1, 2, 4, 5, 6],
  [0, 3, 1, 7, 8, 4],
  [0, 2, 3, 6, 9, 7],
  [1, 3, 2, 8, 9, 5],
];
const HEX_FACES = [
  [0, 1, 2, 3],
  [4, 5, 6, 7],
  [0, 1, 5, 4],
  [1, 2, 6, 5],
  [2, 3, 7, 6],
  [3, 0, 4, 7],
];
const TET_FACES = [
  [0, 1, 2],
  [0, 3, 1],
  [0, 2, 3],
  [1, 3, 2],
];

function fieldValues(mesh: FemResult, result: FemResult | null, field: FieldId): number[] | null {
  const n = mesh.nodeIds?.length ?? 0;
  if (!n) return null;
  if (!result?.ok || result.kind !== "solve") return null;
  const out = new Array<number>(n).fill(0);
  if (field === "vm" && result.vonMises) return result.vonMises;
  if (field === "u" && result.u) {
    for (let i = 0; i < n; i++) {
      const x = result.u[3 * i] ?? 0;
      const y = result.u[3 * i + 1] ?? 0;
      const z = result.u[3 * i + 2] ?? 0;
      out[i] = Math.sqrt(x * x + y * y + z * z);
    }
    return out;
  }
  if ((field === "ux" || field === "uy" || field === "uz") && result.u) {
    const c = field === "ux" ? 0 : field === "uy" ? 1 : 2;
    for (let i = 0; i < n; i++) out[i] = result.u[3 * i + c] ?? 0;
    return out;
  }
  if (result.stress) {
    const map: Record<string, number> = { sxx: 0, syy: 1, szz: 2, sxy: 3 };
    const c = map[field] ?? 0;
    for (let i = 0; i < n; i++) out[i] = result.stress[6 * i + c] ?? 0;
    return out;
  }
  return null;
}

function nodeIndex(mesh: FemResult) {
  const map = new Map<number, number>();
  mesh.nodeIds?.forEach((id, i) => map.set(id, i));
  return map;
}

export function MeshViewer({ mesh, result, field, deformed, scale }: Props) {
  const wrapRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const state = useRef<{
    renderer?: THREE.WebGLRenderer;
    scene?: THREE.Scene;
    camera?: THREE.PerspectiveCamera;
    controls?: OrbitControls;
    solid?: THREE.Mesh;
    edges?: THREE.LineSegments;
    raf?: number;
  }>({});

  useEffect(() => {
    const canvas = canvasRef.current;
    const wrap = wrapRef.current;
    if (!canvas || !wrap) return;

    const renderer = new THREE.WebGLRenderer({
      canvas,
      antialias: true,
      alpha: false,
      powerPreference: "high-performance",
    });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.setClearColor(0x0c1014, 1);
    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(42, 1, 0.05, 5000);
    camera.position.set(80, 50, 110);
    const controls = new OrbitControls(camera, canvas);
    controls.enableDamping = true;
    controls.dampingFactor = 0.08;
    controls.target.set(0, 0, 0);

    scene.add(new THREE.AmbientLight(0xb8c4ce, 0.85));
    const key = new THREE.DirectionalLight(0xf2f4f6, 0.9);
    key.position.set(0.6, 1, 0.4);
    scene.add(key);
    const fill = new THREE.DirectionalLight(0x7f93a3, 0.35);
    fill.position.set(-0.8, 0.2, -0.6);
    scene.add(fill);

    const grid = new THREE.GridHelper(120, 12, 0x2a3340, 0x1a212a);
    grid.position.y = 0;
    scene.add(grid);

    state.current = { renderer, scene, camera, controls };

    const resize = () => {
      const w = wrap.clientWidth || 1;
      const h = wrap.clientHeight || 1;
      renderer.setSize(w, h, false);
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
    };
    resize();
    const ro = new ResizeObserver(resize);
    ro.observe(wrap);

    const loop = () => {
      controls.update();
      renderer.render(scene, camera);
      state.current.raf = requestAnimationFrame(loop);
    };
    loop();

    return () => {
      cancelAnimationFrame(state.current.raf ?? 0);
      ro.disconnect();
      controls.dispose();
      renderer.dispose();
      state.current = {};
    };
  }, []);

  useEffect(() => {
    const { scene, camera, controls } = state.current;
    if (!scene || !camera || !controls) return;

    if (state.current.solid) {
      scene.remove(state.current.solid);
      state.current.solid.geometry.dispose();
      (state.current.solid.material as THREE.Material).dispose();
      state.current.solid = undefined;
    }
    if (state.current.edges) {
      scene.remove(state.current.edges);
      state.current.edges.geometry.dispose();
      (state.current.edges.material as THREE.Material).dispose();
      state.current.edges = undefined;
    }

    if (!mesh?.ok || !mesh.coords || !mesh.nodeIds || !mesh.elements) return;

    const n = mesh.nodeIds.length;
    const idxOf = nodeIndex(mesh);
    const rest: number[] = mesh.coords.slice();
    const pos = rest.slice();
    const u = result?.u;
    if (deformed && u && u.length >= n * 3) {
      for (let i = 0; i < n; i++) {
        pos[3 * i] = rest[3 * i] + scale * u[3 * i];
        pos[3 * i + 1] = rest[3 * i + 1] + scale * u[3 * i + 1];
        pos[3 * i + 2] = rest[3 * i + 2] + scale * u[3 * i + 2];
      }
    }

    const values = fieldValues(mesh, result, field);
    let vmin = 0;
    let vmax = 1;
    if (values && values.length) {
      vmin = Math.min(...values);
      vmax = Math.max(...values);
      if (vmax - vmin < 1e-18) {
        vmax = vmin + 1;
      }
    }

    const positions: number[] = [];
    const colors: number[] = [];
    const pushTri = (a: number, b: number, c: number) => {
      const ids = [a, b, c];
      for (const i of ids) {
        positions.push(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
        if (values) {
          const t = (values[i] - vmin) / (vmax - vmin);
          const [r, g, bl] = sampleColor(t);
          colors.push(r, g, bl);
        } else {
          colors.push(0.55, 0.62, 0.68);
        }
      }
    };

    const nodeColor = (i: number): [number, number, number] => {
      if (values) {
        const t = (values[i] - vmin) / (vmax - vmin);
        return sampleColor(t);
      }
      return [0.72, 0.58, 0.38];
    };

    const mix3 = (
      a: [number, number, number],
      b: [number, number, number],
      c: [number, number, number],
      n1: number,
      n2: number,
      n3: number,
    ): [number, number, number] => [
      a[0] * n1 + b[0] * n2 + c[0] * n3,
      a[1] * n1 + b[1] * n2 + c[1] * n3,
      a[2] * n1 + b[2] * n2 + c[2] * n3,
    ];

    const extraEdges: number[] = [];

    const emitBeam = (
      i0: number,
      i1: number,
      i2: number,
      wa: number,
      hb: number,
      n1h: number[],
      quadratic: boolean,
    ) => {
      const nseg = quadratic ? 10 : 6;
      const p0 = [pos[3 * i0], pos[3 * i0 + 1], pos[3 * i0 + 2]];
      const p1 = [pos[3 * i1], pos[3 * i1 + 1], pos[3 * i1 + 2]];
      const p2 = quadratic
        ? [pos[3 * i2], pos[3 * i2 + 1], pos[3 * i2 + 2]]
        : p0;
      const c0 = nodeColor(i0);
      const c1 = nodeColor(i1);
      const c2 = quadratic ? nodeColor(i2) : c0;
      const hx = n1h[0] ?? 0;
      const hy = n1h[1] ?? 0;
      const hz = n1h[2] ?? -1;
      const ha = Math.max(wa, 1e-6) * 0.5;
      const hbh = Math.max(hb, 1e-6) * 0.5;
      type St = { c: number[]; n1: number[]; n2: number[]; col: [number, number, number] };
      const st: St[] = [];
      for (let s = 0; s <= nseg; s++) {
        const xi = -1 + (2 * s) / nseg;
        let N1: number, N2: number, N3: number, d1: number, d2: number, d3: number;
        if (quadratic) {
          N1 = 0.5 * xi * (xi - 1);
          N2 = 0.5 * xi * (xi + 1);
          N3 = 1 - xi * xi;
          d1 = xi - 0.5;
          d2 = xi + 0.5;
          d3 = -2 * xi;
        } else {
          N1 = 0.5 * (1 - xi);
          N2 = 0.5 * (1 + xi);
          N3 = 0;
          d1 = -0.5;
          d2 = 0.5;
          d3 = 0;
        }
        const cx = N1 * p0[0] + N2 * p1[0] + N3 * p2[0];
        const cy = N1 * p0[1] + N2 * p1[1] + N3 * p2[1];
        const cz = N1 * p0[2] + N2 * p1[2] + N3 * p2[2];
        let tx = d1 * p0[0] + d2 * p1[0] + d3 * p2[0];
        let ty = d1 * p0[1] + d2 * p1[1] + d3 * p2[1];
        let tz = d1 * p0[2] + d2 * p1[2] + d3 * p2[2];
        const tl = Math.hypot(tx, ty, tz) || 1;
        tx /= tl;
        ty /= tl;
        tz /= tl;
        let x1 = hx - (hx * tx + hy * ty + hz * tz) * tx;
        let y1 = hy - (hx * tx + hy * ty + hz * tz) * ty;
        let z1 = hz - (hx * tx + hy * ty + hz * tz) * tz;
        let n1l = Math.hypot(x1, y1, z1);
        if (n1l < 1e-8) {
          const ax = Math.abs(tz) < 0.9 ? 0 : 1;
          const ay = 0;
          const az = Math.abs(tz) < 0.9 ? -1 : 0;
          x1 = ax - (ax * tx + ay * ty + az * tz) * tx;
          y1 = ay - (ax * tx + ay * ty + az * tz) * ty;
          z1 = az - (ax * tx + ay * ty + az * tz) * tz;
          n1l = Math.hypot(x1, y1, z1) || 1;
        }
        x1 /= n1l;
        y1 /= n1l;
        z1 /= n1l;
        const x2 = ty * z1 - tz * y1;
        const y2 = tz * x1 - tx * z1;
        const z2 = tx * y1 - ty * x1;
        st.push({
          c: [cx, cy, cz],
          n1: [x1, y1, z1],
          n2: [x2, y2, z2],
          col: mix3(c0, c1, c2, N1, N2, N3),
        });
      }
      const corner = (s: St, sy: number, sz: number) => [
        s.c[0] + sy * ha * s.n1[0] + sz * hbh * s.n2[0],
        s.c[1] + sy * ha * s.n1[1] + sz * hbh * s.n2[1],
        s.c[2] + sy * ha * s.n1[2] + sz * hbh * s.n2[2],
      ];
      const signs: [number, number][] = [
        [-1, -1],
        [1, -1],
        [1, 1],
        [-1, 1],
      ];
      const pushRaw = (
        pa: number[],
        pb: number[],
        pc: number[],
        ca: [number, number, number],
        cb: [number, number, number],
        cc: [number, number, number],
      ) => {
        positions.push(pa[0], pa[1], pa[2], pb[0], pb[1], pb[2], pc[0], pc[1], pc[2]);
        colors.push(ca[0], ca[1], ca[2], cb[0], cb[1], cb[2], cc[0], cc[1], cc[2]);
      };
      for (let s = 0; s < nseg; s++) {
        const a = st[s];
        const b = st[s + 1];
        extraEdges.push(a.c[0], a.c[1], a.c[2], b.c[0], b.c[1], b.c[2]);
        for (let k = 0; k < 4; k++) {
          const k2 = (k + 1) % 4;
          const a0 = corner(a, signs[k][0], signs[k][1]);
          const a1 = corner(a, signs[k2][0], signs[k2][1]);
          const b0 = corner(b, signs[k][0], signs[k][1]);
          const b1 = corner(b, signs[k2][0], signs[k2][1]);
          pushRaw(a0, b0, b1, a.col, b.col, b.col);
          pushRaw(a0, b1, a1, a.col, b.col, a.col);
        }
      }
      const cap = (s: St, flip: boolean) => {
        const q = signs.map(([sy, sz]) => corner(s, sy, sz));
        if (flip) {
          pushRaw(q[0], q[3], q[2], s.col, s.col, s.col);
          pushRaw(q[0], q[2], q[1], s.col, s.col, s.col);
        } else {
          pushRaw(q[0], q[1], q[2], s.col, s.col, s.col);
          pushRaw(q[0], q[2], q[3], s.col, s.col, s.col);
        }
      };
      cap(st[0], true);
      cap(st[st.length - 1], false);
    };

    const edgeSet = new Set<string>();
    const addEdge = (a: number, b: number) => {
      const lo = Math.min(a, b);
      const hi = Math.max(a, b);
      edgeSet.add(`${lo}-${hi}`);
    };

    for (const el of mesh.elements) {
      const loc = el.nodes.map((id) => idxOf.get(id)).filter((v): v is number => v !== undefined);
      const t = el.type;
      if ((t === "C3D8" || t === "C3D8R" || t === "C3D8I") && loc.length >= 8) {
        for (const f of HEX_FACES) {
          pushTri(loc[f[0]], loc[f[1]], loc[f[2]]);
          pushTri(loc[f[0]], loc[f[2]], loc[f[3]]);
          addEdge(loc[f[0]], loc[f[1]]);
          addEdge(loc[f[1]], loc[f[2]]);
          addEdge(loc[f[2]], loc[f[3]]);
          addEdge(loc[f[3]], loc[f[0]]);
        }
      } else if (t === "C3D4" && loc.length >= 4) {
        for (const f of TET_FACES) {
          pushTri(loc[f[0]], loc[f[1]], loc[f[2]]);
          addEdge(loc[f[0]], loc[f[1]]);
          addEdge(loc[f[1]], loc[f[2]]);
          addEdge(loc[f[2]], loc[f[0]]);
        }
      } else if ((t === "C3D6" || t === "C3D15") && loc.length >= 6) {
        const bot = [loc[0], loc[1], loc[2]];
        const top = [loc[3], loc[4], loc[5]];
        pushTri(bot[0], bot[1], bot[2]);
        pushTri(top[0], top[2], top[1]);
        addEdge(bot[0], bot[1]);
        addEdge(bot[1], bot[2]);
        addEdge(bot[2], bot[0]);
        addEdge(top[0], top[1]);
        addEdge(top[1], top[2]);
        addEdge(top[2], top[0]);
        addEdge(bot[0], top[0]);
        addEdge(bot[1], top[1]);
        addEdge(bot[2], top[2]);
        if (t === "C3D15" && loc.length >= 15) {
          addEdge(bot[0], loc[6]);
          addEdge(loc[6], bot[1]);
          addEdge(top[0], loc[9]);
          addEdge(loc[9], top[1]);
        }
      } else if (
        (t === "CPS4" ||
          t === "CPE4" ||
          t === "S4" ||
          t === "S4R" ||
          t === "CAX4" ||
          t === "CAX4R" ||
          t === "M3D4" ||
          t === "M3D4R") &&
        loc.length >= 4
      ) {
        pushTri(loc[0], loc[1], loc[2]);
        pushTri(loc[0], loc[2], loc[3]);
        addEdge(loc[0], loc[1]);
        addEdge(loc[1], loc[2]);
        addEdge(loc[2], loc[3]);
        addEdge(loc[3], loc[0]);
      } else if ((t === "CPS3" || t === "CPE3" || t === "S3" || t === "S3R") && loc.length >= 3) {
        pushTri(loc[0], loc[1], loc[2]);
        addEdge(loc[0], loc[1]);
        addEdge(loc[1], loc[2]);
        addEdge(loc[2], loc[0]);
      } else if ((t === "C3D20" || t === "C3D20R") && loc.length >= 20) {
        for (const f of HEX20_FACES) {
          const [c0, c1, c2, c3, m01, m12, m23, m30] = f.map((i) => loc[i]);
          pushTri(c0, m01, m30);
          pushTri(m01, c1, m12);
          pushTri(m12, c2, m23);
          pushTri(m23, c3, m30);
          pushTri(m01, m12, m23);
          pushTri(m01, m23, m30);
          addEdge(c0, m01);
          addEdge(m01, c1);
          addEdge(c1, m12);
          addEdge(m12, c2);
          addEdge(c2, m23);
          addEdge(m23, c3);
          addEdge(c3, m30);
          addEdge(m30, c0);
        }
      } else if (t === "C3D10" && loc.length >= 10) {
        for (const f of TET10_FACES) {
          const [c0, c1, c2, m01, m12, m20] = f.map((i) => loc[i]);
          pushTri(c0, m01, m20);
          pushTri(m01, c1, m12);
          pushTri(m20, m12, c2);
          pushTri(m01, m12, m20);
          addEdge(c0, m01);
          addEdge(m01, c1);
          addEdge(c1, m12);
          addEdge(m12, c2);
          addEdge(c2, m20);
          addEdge(m20, c0);
        }
      } else if (
        (t === "CPS8" ||
          t === "CPE8" ||
          t === "CPS8R" ||
          t === "CPE8R" ||
          t === "S8" ||
          t === "S8R" ||
          t === "CAX8" ||
          t === "CAX8R" ||
          t === "M3D8") &&
        loc.length >= 8
      ) {
        const [c0, c1, c2, c3, m01, m12, m23, m30] = loc;
        pushTri(c0, m01, m30);
        pushTri(m01, c1, m12);
        pushTri(m12, c2, m23);
        pushTri(m23, c3, m30);
        pushTri(m01, m12, m23);
        pushTri(m01, m23, m30);
        addEdge(c0, m01);
        addEdge(m01, c1);
        addEdge(c1, m12);
        addEdge(m12, c2);
        addEdge(c2, m23);
        addEdge(m23, c3);
        addEdge(c3, m30);
        addEdge(m30, c0);
      } else if ((t === "CPS6" || t === "CPE6" || t === "S6") && loc.length >= 6) {
        const [c0, c1, c2, m01, m12, m20] = loc;
        pushTri(c0, m01, m20);
        pushTri(m01, c1, m12);
        pushTri(m20, m12, c2);
        pushTri(m01, m12, m20);
        addEdge(c0, m01);
        addEdge(m01, c1);
        addEdge(c1, m12);
        addEdge(m12, c2);
        addEdge(c2, m20);
        addEdge(m20, c0);
      } else if ((t === "B32" || t === "B32R") && loc.length >= 3) {
        emitBeam(loc[0], loc[1], loc[2], el.secA ?? 8, el.secB ?? 8, el.n1 ?? [0, 0, -1], true);
      } else if ((t === "B31" || t === "B31R") && loc.length >= 2) {
        emitBeam(loc[0], loc[1], loc[0], el.secA ?? 8, el.secB ?? 8, el.n1 ?? [0, 0, -1], false);
      }
    }

    if (positions.length === 0) return;

    const geo = new THREE.BufferGeometry();
    geo.setAttribute("position", new THREE.Float32BufferAttribute(positions, 3));
    geo.setAttribute("color", new THREE.Float32BufferAttribute(colors, 3));
    geo.computeVertexNormals();
    const mat = new THREE.MeshLambertMaterial({
      vertexColors: true,
      side: THREE.DoubleSide,
      polygonOffset: true,
      polygonOffsetFactor: 1,
      polygonOffsetUnits: 1,
    });
    const solid = new THREE.Mesh(geo, mat);
    scene.add(solid);
    state.current.solid = solid;

    const epos: number[] = [...extraEdges];
    for (const key of edgeSet) {
      const [a, b] = key.split("-").map(Number);
      epos.push(pos[3 * a], pos[3 * a + 1], pos[3 * a + 2], pos[3 * b], pos[3 * b + 1], pos[3 * b + 2]);
    }
    const ego = new THREE.BufferGeometry();
    ego.setAttribute("position", new THREE.Float32BufferAttribute(epos, 3));
    const edges = new THREE.LineSegments(
      ego,
      new THREE.LineBasicMaterial({ color: 0x0a0c0e, transparent: true, opacity: 0.55 }),
    );
    scene.add(edges);
    state.current.edges = edges;

    const lastFit = state.current as typeof state.current & { fitKey?: string };
    geo.computeBoundingBox();
    const bb = geo.boundingBox;
    const fitKey = `${mesh.nnode}-${mesh.nelem}-${mesh.heading}`;
    if (bb && lastFit.fitKey !== fitKey) {
      lastFit.fitKey = fitKey;
      const size = new THREE.Vector3();
      bb.getSize(size);
      const center = new THREE.Vector3();
      bb.getCenter(center);
      const maxDim = Math.max(size.x, size.y, size.z, 1);
      const dist = maxDim * 1.85;
      camera.near = maxDim / 200;
      camera.far = maxDim * 40;
      camera.updateProjectionMatrix();
      const is2d = (mesh.dim ?? 3) === 2 || size.z < maxDim * 0.02;
      if (is2d) {
        camera.position.set(center.x, center.y, center.z + dist);
        camera.up.set(0, 1, 0);
      } else {
        camera.position.set(center.x + dist * 0.7, center.y + dist * 0.45, center.z + dist * 0.7);
      }
      controls.target.copy(center);
      controls.update();
    }
  }, [mesh, result, field, deformed, scale]);

  return (
    <div ref={wrapRef} className="relative h-full min-h-[240px] w-full overflow-hidden bg-viewport">
      <canvas ref={canvasRef} className="block h-full w-full touch-none" />
      {!mesh?.ok && (
        <div className="pointer-events-none absolute inset-0 flex items-center justify-center">
          <p className="text-sm text-muted">Netz wird gelesen…</p>
        </div>
      )}
    </div>
  );
}

export function fieldRange(
  mesh: FemResult | null,
  result: FemResult | null,
  field: FieldId,
): { min: number; max: number } | null {
  if (!mesh) return null;
  const v = fieldValues(mesh, result, field);
  if (!v || !v.length) return null;
  return { min: Math.min(...v), max: Math.max(...v) };
}
