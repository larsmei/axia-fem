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

    const edgeSet = new Set<string>();
    const addEdge = (a: number, b: number) => {
      const lo = Math.min(a, b);
      const hi = Math.max(a, b);
      edgeSet.add(`${lo}-${hi}`);
    };

    for (const el of mesh.elements) {
      const loc = el.nodes.map((id) => idxOf.get(id)).filter((v): v is number => v !== undefined);
      const t = el.type;
      if ((t === "C3D8" || t === "C3D8R") && loc.length >= 8) {
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
      } else if ((t === "CPS4" || t === "CPE4") && loc.length >= 4) {
        pushTri(loc[0], loc[1], loc[2]);
        pushTri(loc[0], loc[2], loc[3]);
        addEdge(loc[0], loc[1]);
        addEdge(loc[1], loc[2]);
        addEdge(loc[2], loc[3]);
        addEdge(loc[3], loc[0]);
      } else if ((t === "CPS3" || t === "CPE3") && loc.length >= 3) {
        pushTri(loc[0], loc[1], loc[2]);
        addEdge(loc[0], loc[1]);
        addEdge(loc[1], loc[2]);
        addEdge(loc[2], loc[0]);
      }
    }

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

    const epos: number[] = [];
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
