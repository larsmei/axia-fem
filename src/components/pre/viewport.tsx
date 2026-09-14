import { useCallback, useEffect, useRef } from "react";
import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { resolveNodes } from "@/lib/pre/export-inp";
import { dist, modelBounds3, nearestNode3, shapeFaces, snapPt } from "@/lib/pre/geometry";
import { usePre } from "@/lib/pre/store";
import type { FaceName, Mesh, MeshNode, Shape, Vec2, Vec3 } from "@/lib/pre/types";

const BG = 0xf4f5f7;
const STEEL = 0x8b929c;
const STEEL_SEL = 0x3f4450;
const WIRE = 0x2a2a32;
const FACE_HIT = 0x18181c;
const DRAFT = 0x18181c;

const HEX_FACES = [
  [0, 1, 2, 3],
  [4, 5, 6, 7],
  [0, 1, 5, 4],
  [1, 2, 6, 5],
  [2, 3, 7, 6],
  [3, 0, 4, 7],
];
const WEDGE_FACES = [
  [0, 1, 2],
  [3, 4, 5],
  [0, 1, 4, 3],
  [1, 2, 5, 4],
  [2, 0, 3, 5],
];
const TET_FACES = [
  [0, 1, 2],
  [0, 3, 1],
  [0, 2, 3],
  [1, 3, 2],
];

function to3(p: { x: number; y: number; z?: number }): THREE.Vector3 {
  return new THREE.Vector3(p.x, p.z ?? 0, p.y);
}

function from3(v: THREE.Vector3): Vec3 {
  return { x: v.x, y: v.z, z: v.y };
}

function disposeObj(o: THREE.Object3D) {
  o.traverse((c) => {
    if (c instanceof THREE.Mesh || c instanceof THREE.LineSegments || c instanceof THREE.Line || c instanceof THREE.Points) {
      c.geometry.dispose();
      const m = c.material;
      if (Array.isArray(m)) m.forEach((x) => x.dispose());
      else m.dispose();
    }
  });
}

function makeQuadGeo(corners: Vec3[]): THREE.BufferGeometry {
  const a = to3(corners[0]),
    b = to3(corners[1]),
    c = to3(corners[2]),
    d = to3(corners[3]);
  const g = new THREE.BufferGeometry();
  g.setAttribute(
    "position",
    new THREE.Float32BufferAttribute(
      [a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z, a.x, a.y, a.z, c.x, c.y, c.z, d.x, d.y, d.z],
      3,
    ),
  );
  g.computeVertexNormals();
  return g;
}

function meshToGeometry(mesh: Mesh, filter?: (n: MeshNode) => boolean): { solid: THREE.BufferGeometry; edges: THREE.BufferGeometry } {
  const byId = new Map(mesh.nodes.map((n) => [n.id, n]));
  const pos: number[] = [];
  const edgeSet = new Set<string>();
  const pushTri = (a: MeshNode, b: MeshNode, c: MeshNode) => {
    const pa = to3(a),
      pb = to3(b),
      pc = to3(c);
    pos.push(pa.x, pa.y, pa.z, pb.x, pb.y, pb.z, pc.x, pc.y, pc.z);
  };
  const addEdge = (a: MeshNode, b: MeshNode) => {
    const lo = Math.min(a.id, b.id);
    const hi = Math.max(a.id, b.id);
    edgeSet.add(`${lo}-${hi}`);
  };
  const emitFace = (ids: number[]) => {
    const pts = ids.map((id) => byId.get(id)).filter((n): n is MeshNode => !!n);
    if (filter && pts.some((p) => !filter(p))) return;
    if (pts.length < 3) return;
    if (pts.length === 3) {
      pushTri(pts[0], pts[1], pts[2]);
      addEdge(pts[0], pts[1]);
      addEdge(pts[1], pts[2]);
      addEdge(pts[2], pts[0]);
      return;
    }
    pushTri(pts[0], pts[1], pts[2]);
    pushTri(pts[0], pts[2], pts[3]);
    for (let i = 0; i < 4; i++) addEdge(pts[i], pts[(i + 1) % 4]);
  };
  for (const el of mesh.elements) {
    if ((el.type === "C3D8" || el.type === "CPS4") && el.nodes.length >= (el.type === "C3D8" ? 8 : 4)) {
      if (el.type === "C3D8") for (const f of HEX_FACES) emitFace(f.map((i) => el.nodes[i]));
      else emitFace(el.nodes.slice(0, 4));
    } else if (el.type === "C3D6" && el.nodes.length >= 6) {
      for (const f of WEDGE_FACES) emitFace(f.map((i) => el.nodes[i]));
    } else if (el.type === "C3D4" && el.nodes.length >= 4) {
      for (const f of TET_FACES) emitFace(f.map((i) => el.nodes[i]));
    } else if (el.nodes.length >= 3) {
      emitFace(el.nodes.slice(0, 3));
    }
  }
  const solid = new THREE.BufferGeometry();
  solid.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
  solid.computeVertexNormals();
  const epos: number[] = [];
  for (const key of edgeSet) {
    const [a, b] = key.split("-").map(Number);
    const na = byId.get(a),
      nb = byId.get(b);
    if (!na || !nb) continue;
    const pa = to3(na),
      pb = to3(nb);
    epos.push(pa.x, pa.y, pa.z, pb.x, pb.y, pb.z);
  }
  const edges = new THREE.BufferGeometry();
  edges.setAttribute("position", new THREE.Float32BufferAttribute(epos, 3));
  return { solid, edges };
}

function solidGeometry(shape: Shape): THREE.BufferGeometry {
  if (shape.kind === "box") {
    const g = new THREE.BoxGeometry(shape.w, shape.d, shape.h);
    g.translate(shape.x + shape.w / 2, shape.z + shape.d / 2, shape.y + shape.h / 2);
    return g;
  }
  if (shape.kind === "rect") {
    const d = Math.max(shape.depth, 0.4);
    const g = new THREE.BoxGeometry(shape.w, d, shape.h);
    g.translate(shape.x + shape.w / 2, d / 2, shape.y + shape.h / 2);
    return g;
  }
  if (shape.kind === "sphere") {
    const g = new THREE.SphereGeometry(shape.r, 32, 20);
    g.translate(shape.cx, shape.cz, shape.cy);
    return g;
  }
  if (shape.kind === "cylinder") {
    const g = new THREE.CylinderGeometry(shape.r, shape.r, shape.height, 32);
    if (shape.axis === "z") {
      /* three Y is FEM Z */
    } else if (shape.axis === "y") {
      g.rotateX(Math.PI / 2);
    } else {
      g.rotateZ(Math.PI / 2);
    }
    g.translate(shape.cx, shape.cz, shape.cy);
    return g;
  }
  if (shape.kind === "circle") {
    const d = Math.max(shape.depth, 0.4);
    const g = new THREE.CylinderGeometry(shape.r, shape.r, d, 32);
    g.translate(shape.cx, d / 2, shape.cy);
    return g;
  }
  const sh = new THREE.Shape();
  shape.points.forEach((p, i) => {
    if (i === 0) sh.moveTo(p.x, p.y);
    else sh.lineTo(p.x, p.y);
  });
  sh.closePath();
  for (const h of shape.holes ?? []) {
    const hole = new THREE.Path();
    hole.absarc(h.cx, h.cy, h.r, 0, Math.PI * 2, true);
    sh.holes.push(hole);
  }
  const g = new THREE.ExtrudeGeometry(sh, { depth: Math.max(shape.depth, 0.4), bevelEnabled: false });
  g.rotateX(-Math.PI / 2);
  return g;
}

export function PreViewport() {
  const wrapRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const st = useRef<{
    renderer?: THREE.WebGLRenderer;
    scene?: THREE.Scene;
    camera?: THREE.PerspectiveCamera;
    controls?: OrbitControls;
    world?: THREE.Group;
    picks?: THREE.Group;
    draft?: THREE.Object3D;
    axesScene?: THREE.Scene;
    axesCam?: THREE.PerspectiveCamera;
    ray?: THREE.Raycaster;
    plane?: THREE.Plane;
    raf?: number;
  }>({});
  const shapes = usePre((s) => s.shapes);
  const mesh = usePre((s) => s.mesh);
  const tool = usePre((s) => s.tool);
  const dim = usePre((s) => s.dim);
  const selectedShapeId = usePre((s) => s.selectedShapeId);
  const selectedFace = usePre((s) => s.selectedFace);
  const selectedNodeIds = usePre((s) => s.selectedNodeIds);
  const restraints = usePre((s) => s.restraints);
  const loads = usePre((s) => s.loads);
  const draft = usePre((s) => s.draft);

  const fit = useCallback((mode: "iso" | "top" | "keep" = "iso") => {
    const { camera, controls } = st.current;
    if (!camera || !controls) return;
    const b = modelBounds3(usePre.getState().shapes);
    const cx = b.x + b.w / 2;
    const cy = b.y + b.h / 2;
    const cz = b.z + b.d / 2;
    const span = Math.max(b.w, b.h, b.d, 30);
    camera.near = Math.max(0.05, span / 400);
    camera.far = span * 50;
    camera.updateProjectionMatrix();
    const t = to3({ x: cx, y: cy, z: cz });
    controls.target.copy(t);
    if (mode === "top") {
      camera.up.set(0, 0, -1);
      camera.position.set(t.x, t.y + span * 1.8, t.z);
    } else {
      camera.up.set(0, 1, 0);
      camera.position.set(t.x + span * 0.95, t.y + span * 0.75, t.z + span * 1.15);
    }
    controls.update();
  }, []);

  useEffect(() => {
    const canvas = canvasRef.current;
    const wrap = wrapRef.current;
    if (!canvas || !wrap) return;
    const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: false, powerPreference: "high-performance" });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.setClearColor(BG, 1);
    renderer.autoClear = false;
    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(42, 1, 0.05, 8000);
    camera.up.set(0, 1, 0);
    const controls = new OrbitControls(camera, canvas);
    controls.enableDamping = true;
    controls.dampingFactor = 0.08;
    controls.mouseButtons = {
      LEFT: THREE.MOUSE.ROTATE,
      MIDDLE: THREE.MOUSE.PAN,
      RIGHT: THREE.MOUSE.ROTATE,
    };
    controls.touches = { ONE: THREE.TOUCH.ROTATE, TWO: THREE.TOUCH.DOLLY_PAN };
    scene.add(new THREE.AmbientLight(0xffffff, 0.92));
    const key = new THREE.DirectionalLight(0xffffff, 0.5);
    key.position.set(0.4, 1, 0.35);
    scene.add(key);
    const fill = new THREE.DirectionalLight(0xffffff, 0.22);
    fill.position.set(-0.6, 0.2, -0.5);
    scene.add(fill);

    const grid = new THREE.GridHelper(200, 20, 0xc5cad1, 0xdce0e5);
    scene.add(grid);
    const axes = new THREE.AxesHelper(24);
    axes.setColors(new THREE.Color(0x18181c), new THREE.Color(0x52525b), new THREE.Color(0x71717a));
    scene.add(axes);

    const world = new THREE.Group();
    const picks = new THREE.Group();
    scene.add(world);
    scene.add(picks);

    const axesScene = new THREE.Scene();
    const triad = new THREE.Group();
    const axisLine = (to: [number, number, number]) => {
      const g = new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(0, 0, 0), new THREE.Vector3(...to)]);
      return new THREE.Line(g, new THREE.LineBasicMaterial({ color: 0x1a1a1a }));
    };
    triad.add(axisLine([1, 0, 0]), axisLine([0, 1, 0]), axisLine([0, 0, 1]));
    const makeLbl = (t: string, p: [number, number, number]) => {
      const c = document.createElement("canvas");
      c.width = 64;
      c.height = 64;
      const ctx = c.getContext("2d")!;
      ctx.fillStyle = "#18181c";
      ctx.font = "28px IBM Plex Sans, sans-serif";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText(t, 32, 32);
      const tex = new THREE.CanvasTexture(c);
      const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, depthTest: false }));
      sp.position.set(...p);
      sp.scale.setScalar(0.45);
      return sp;
    };
    triad.add(makeLbl("X", [1.25, 0, 0]), makeLbl("Z", [0, 1.25, 0]), makeLbl("Y", [0, 0, 1.25]));
    axesScene.add(triad);
    const axesCam = new THREE.PerspectiveCamera(50, 1, 0.1, 10);

    st.current = {
      renderer,
      scene,
      camera,
      controls,
      world,
      picks,
      axesScene,
      axesCam,
      ray: new THREE.Raycaster(),
      plane: new THREE.Plane(new THREE.Vector3(0, 1, 0), 0),
    };

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
      const w = wrap.clientWidth || 1;
      const h = wrap.clientHeight || 1;
      renderer.setViewport(0, 0, w, h);
      renderer.setScissorTest(false);
      renderer.clear();
      renderer.render(scene, camera);
      axesCam.up.copy(camera.up);
      axesCam.position.copy(camera.position).sub(controls.target);
      if (axesCam.position.lengthSq() < 1e-12) axesCam.position.set(0.7, 0.45, 0.7);
      axesCam.position.setLength(2.4);
      axesCam.lookAt(0, 0, 0);
      const aw = 72,
        ah = 72;
      renderer.clearDepth();
      renderer.setScissorTest(true);
      renderer.setScissor(w - aw - 10, 10, aw, ah);
      renderer.setViewport(w - aw - 10, 10, aw, ah);
      renderer.render(axesScene, axesCam);
      renderer.setScissorTest(false);
      renderer.setViewport(0, 0, w, h);
      st.current.raf = requestAnimationFrame(loop);
    };
    loop();
    fit(usePre.getState().dim === "2d" ? "top" : "iso");

    return () => {
      cancelAnimationFrame(st.current.raf ?? 0);
      ro.disconnect();
      controls.dispose();
      renderer.dispose();
      disposeObj(world);
      disposeObj(picks);
      st.current = {};
    };
  }, [fit]);

  useEffect(() => {
    const { world, picks } = st.current;
    if (!world || !picks) return;
    while (world.children.length) {
      const c = world.children[0];
      world.remove(c);
      disposeObj(c);
    }
    while (picks.children.length) {
      const c = picks.children[0];
      picks.remove(c);
      disposeObj(c);
    }

    const state = usePre.getState();
    const matFor = (sel: boolean, opacity = 0.92) =>
      new THREE.MeshLambertMaterial({
        color: sel ? STEEL_SEL : STEEL,
        transparent: opacity < 0.99,
        opacity,
        side: THREE.DoubleSide,
      });

    if (state.mesh && state.mesh.elements.length) {
      const { solid, edges } = meshToGeometry(state.mesh);
      const meshObj = new THREE.Mesh(solid, matFor(false, 0.96));
      const wire = new THREE.LineSegments(edges, new THREE.LineBasicMaterial({ color: WIRE, transparent: true, opacity: 0.45 }));
      world.add(meshObj, wire);
      if (state.selectedNodeIds.length) {
        const sel = new THREE.BufferGeometry();
        const pts: number[] = [];
        const set = new Set(state.selectedNodeIds);
        for (const n of state.mesh.nodes) {
          if (!set.has(n.id)) continue;
          const p = to3(n);
          pts.push(p.x, p.y, p.z);
        }
        sel.setAttribute("position", new THREE.Float32BufferAttribute(pts, 3));
        world.add(new THREE.Points(sel, new THREE.PointsMaterial({ color: 0x09090b, size: 7, sizeAttenuation: false })));
      }
    } else {
      for (const sh of state.shapes) {
        const geo = solidGeometry(sh);
        const m = new THREE.Mesh(geo, matFor(sh.id === state.selectedShapeId, 0.88));
        world.add(m);
        world.add(new THREE.LineSegments(new THREE.EdgesGeometry(geo), new THREE.LineBasicMaterial({ color: WIRE, opacity: 0.55, transparent: true })));
      }
    }

    for (const sh of state.shapes) {
      for (const f of shapeFaces(sh)) {
        const g = makeQuadGeo(f.corners);
        const m = new THREE.Mesh(g, new THREE.MeshBasicMaterial({ side: THREE.DoubleSide, transparent: true, opacity: 0 }));
        m.userData = { shapeId: sh.id, face: f.face };
        picks.add(m);
        if (state.selectedFace?.shapeId === sh.id && state.selectedFace.face === f.face) {
          const hl = new THREE.Mesh(
            g.clone(),
            new THREE.MeshBasicMaterial({ color: FACE_HIT, transparent: true, opacity: 0.22, side: THREE.DoubleSide, depthWrite: false }),
          );
          world.add(hl);
        }
      }
      if (sh.kind === "cylinder" || sh.kind === "sphere") {
        const geo = solidGeometry(sh);
        const m = new THREE.Mesh(geo, new THREE.MeshBasicMaterial({ visible: false, side: THREE.DoubleSide }));
        m.userData = { shapeId: sh.id, face: "lateral" };
        picks.add(m);
      }
    }

    if (state.mesh) {
      for (const ld of state.loads) {
        if (ld.kind !== "force") continue;
        const ids = resolveNodes(state.mesh, state.shapes, ld.target);
        if (!ids.length) continue;
        let cx = 0,
          cy = 0,
          cz = 0;
        let n = 0;
        for (const id of ids) {
          const nd = state.mesh.nodes.find((x) => x.id === id);
          if (!nd) continue;
          cx += nd.x;
          cy += nd.y;
          cz += nd.z ?? 0;
          n++;
        }
        if (!n) continue;
        cx /= n;
        cy /= n;
        cz /= n;
        const mag = Math.hypot(ld.fx, ld.fy, ld.fz) || 1;
        const dir = to3({ x: ld.fx / mag, y: ld.fy / mag, z: ld.fz / mag });
        if (dir.lengthSq() < 1e-12) continue;
        const origin = to3({ x: cx, y: cy, z: cz });
        world.add(new THREE.ArrowHelper(dir.normalize(), origin, 16, 0x18181c, 5, 3));
      }
      const pinIds: MeshNode[] = [];
      const seen = new Set<number>();
      for (const r of state.restraints) {
        const ids = resolveNodes(state.mesh, state.shapes, r.target);
        const step = Math.max(1, Math.floor(ids.length / 48));
        for (let i = 0; i < ids.length; i += step) {
          if (seen.has(ids[i])) continue;
          seen.add(ids[i]);
          const nd = state.mesh.nodes.find((x) => x.id === ids[i]);
          if (nd) pinIds.push(nd);
        }
      }
      if (pinIds.length) {
        const pin = new THREE.InstancedMesh(
          new THREE.ConeGeometry(1.6, 4.2, 6),
          new THREE.MeshLambertMaterial({ color: 0x18181c }),
          pinIds.length,
        );
        const dummy = new THREE.Object3D();
        pinIds.forEach((nd, i) => {
          const p = to3(nd);
          dummy.position.set(p.x, p.y - 2, p.z);
          dummy.rotation.set(Math.PI, 0, 0);
          dummy.updateMatrix();
          pin.setMatrixAt(i, dummy.matrix);
        });
        pin.instanceMatrix.needsUpdate = true;
        world.add(pin);
      }
    }

    const d = state.draft;
    if (d?.tool === "box" || d?.tool === "rect") {
      const x0 = Math.min(d.a.x, d.b.x),
        y0 = Math.min(d.a.y, d.b.y);
      const w = Math.abs(d.b.x - d.a.x) || 0.01,
        h = Math.abs(d.b.y - d.a.y) || 0.01;
      const depth = state.dim === "3d" ? state.defaultDepth : 0.4;
      const g = new THREE.BoxGeometry(w, depth, h);
      g.translate(x0 + w / 2, depth / 2, y0 + h / 2);
      world.add(
        new THREE.Mesh(g, new THREE.MeshBasicMaterial({ color: DRAFT, transparent: true, opacity: 0.18, depthWrite: false })),
      );
      world.add(new THREE.LineSegments(new THREE.EdgesGeometry(g), new THREE.LineBasicMaterial({ color: DRAFT })));
    }
    if (d?.tool === "cylinder" || d?.tool === "circle" || d?.tool === "sphere" || d?.tool === "hole") {
      const r = Math.max(d.r, 0.01);
      const depth = d.tool === "sphere" ? r * 2 : state.dim === "3d" ? state.defaultDepth : 0.4;
      const g =
        d.tool === "sphere"
          ? new THREE.SphereGeometry(r, 24, 16)
          : new THREE.CylinderGeometry(r, r, depth, 28);
      if (d.tool === "sphere") g.translate(d.c.x, r, d.c.y);
      else g.translate(d.c.x, depth / 2, d.c.y);
      world.add(new THREE.Mesh(g, new THREE.MeshBasicMaterial({ color: DRAFT, transparent: true, opacity: 0.16, depthWrite: false })));
      world.add(new THREE.LineSegments(new THREE.EdgesGeometry(g), new THREE.LineBasicMaterial({ color: DRAFT })));
    }
    if (d?.tool === "polygon" && d.points.length) {
      const pts = d.points.map((p) => to3({ x: p.x, y: p.y, z: 0 }));
      const g = new THREE.BufferGeometry().setFromPoints(pts);
      world.add(new THREE.Line(g, new THREE.LineBasicMaterial({ color: DRAFT })));
    }
  }, [shapes, mesh, selectedShapeId, selectedFace, selectedNodeIds, restraints, loads, draft, dim]);

  useEffect(() => {
    const { controls } = st.current;
    if (!controls) return;
    const draw = tool === "rect" || tool === "circle" || tool === "polygon" || tool === "hole" || tool === "box" || tool === "cylinder" || tool === "sphere";
    controls.mouseButtons.LEFT = draw ? (-1 as THREE.MOUSE) : THREE.MOUSE.ROTATE;
  }, [tool]);

  useEffect(() => {
    fit(dim === "2d" ? "top" : "iso");
  }, [dim, shapes.length, fit]);

  const hitPlane = (e: React.PointerEvent): Vec2 | null => {
    const { camera, ray, plane, renderer } = st.current;
    const canvas = canvasRef.current;
    if (!camera || !ray || !plane || !canvas || !renderer) return null;
    const r = canvas.getBoundingClientRect();
    const ndc = new THREE.Vector2(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
    ray.setFromCamera(ndc, camera);
    const out = new THREE.Vector3();
    if (!ray.ray.intersectPlane(plane, out)) return null;
    const p = from3(out);
    return snapPt({ x: p.x, y: p.y }, 1);
  };

  const hitPick = (e: React.PointerEvent): { shapeId: string; face?: FaceName } | null => {
    const { camera, ray, picks } = st.current;
    const canvas = canvasRef.current;
    if (!camera || !ray || !picks || !canvas) return null;
    const r = canvas.getBoundingClientRect();
    const ndc = new THREE.Vector2(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
    ray.setFromCamera(ndc, camera);
    const hits = ray.intersectObjects(picks.children, false);
    const h = hits[0];
    if (!h?.object.userData?.shapeId) return null;
    return { shapeId: h.object.userData.shapeId, face: h.object.userData.face };
  };

  const onDown = (e: React.PointerEvent) => {
    if (e.button === 1 || e.button === 2) return;
    const pre = usePre.getState();
    const wp = hitPlane(e);
    if (pre.tool === "box" || pre.tool === "rect") {
      if (!wp) return;
      (e.target as HTMLElement).setPointerCapture(e.pointerId);
      pre.setDraft({ tool: pre.tool, a: wp, b: wp });
      return;
    }
    if (pre.tool === "cylinder" || pre.tool === "circle" || pre.tool === "sphere" || pre.tool === "hole") {
      if (!wp) return;
      (e.target as HTMLElement).setPointerCapture(e.pointerId);
      pre.setDraft({ tool: pre.tool, c: wp, r: 0 });
      return;
    }
    if (pre.tool === "polygon") {
      if (!wp) return;
      const d = pre.draft?.tool === "polygon" ? pre.draft.points : [];
      if (d.length >= 3 && dist(wp, d[0]) < 3) {
        pre.addPolygon(d);
        pre.setDraft(null);
        pre.setTool("select");
        return;
      }
      pre.setDraft({ tool: "polygon", points: [...d, wp] });
      return;
    }
    if (pre.tool === "node" && pre.mesh && st.current.camera && st.current.ray) {
      const canvas = canvasRef.current!;
      const r = canvas.getBoundingClientRect();
      const ndc = new THREE.Vector2(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
      st.current.ray.setFromCamera(ndc, st.current.camera);
      const origin = from3(st.current.ray.ray.origin);
      const dirV = from3(st.current.ray.ray.direction);
      const distCam = st.current.camera.position.distanceTo(st.current.controls?.target ?? new THREE.Vector3());
      const n = nearestNode3(pre.mesh, origin, dirV, Math.max(1.5, distCam * 0.02));
      if (n) pre.toggleNode(n.id, e.shiftKey);
      return;
    }
    const hit = hitPick(e);
    if (hit) {
      if (pre.tool === "face" && hit.face) {
        pre.selectFace({ shapeId: hit.shapeId, face: hit.face });
      } else {
        pre.selectShape(hit.shapeId);
        if (hit.face) pre.selectFace({ shapeId: hit.shapeId, face: hit.face });
      }
    } else {
      pre.selectShape(null);
      pre.selectFace(null);
    }
  };

  const onMove = (e: React.PointerEvent) => {
    const pre = usePre.getState();
    const wp = hitPlane(e);
    if (!wp || !pre.draft) return;
    if (pre.draft.tool === "rect" || pre.draft.tool === "box") pre.setDraft({ tool: pre.draft.tool, a: pre.draft.a, b: wp });
    if (pre.draft.tool === "circle" || pre.draft.tool === "cylinder" || pre.draft.tool === "sphere" || pre.draft.tool === "hole") {
      pre.setDraft({ tool: pre.draft.tool, c: pre.draft.c, r: dist(pre.draft.c, wp) });
    }
  };

  const onUp = () => {
    const pre = usePre.getState();
    const d = pre.draft;
    if (d?.tool === "rect" || d?.tool === "box") {
      const x = Math.min(d.a.x, d.b.x);
      const y = Math.min(d.a.y, d.b.y);
      const w = Math.abs(d.b.x - d.a.x);
      const h = Math.abs(d.b.y - d.a.y);
      if (d.tool === "box" || pre.dim === "3d") pre.addBox(x, y, 0, w, h, pre.defaultDepth);
      else pre.addRect(x, y, w, h);
      pre.setDraft(null);
      pre.setTool("select");
    }
    if (d?.tool === "circle" || d?.tool === "cylinder") {
      if (d.tool === "cylinder" || pre.dim === "3d") pre.addCylinder(d.c.x, d.c.y, pre.defaultDepth / 2, d.r, pre.defaultDepth, "z");
      else pre.addCircle(d.c.x, d.c.y, d.r);
      pre.setDraft(null);
      pre.setTool("select");
    }
    if (d?.tool === "sphere") {
      pre.addSphere(d.c.x, d.c.y, d.r, d.r);
      pre.setDraft(null);
      pre.setTool("select");
    }
    if (d?.tool === "hole") {
      pre.addHole(d.c.x, d.c.y, d.r);
      pre.setDraft(null);
      pre.setTool("select");
    }
  };

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      const t = e.target as HTMLElement | null;
      if (t && (t.tagName === "INPUT" || t.tagName === "TEXTAREA" || t.tagName === "SELECT")) return;
      const pre = usePre.getState();
      if (e.key === "Escape") {
        pre.setDraft(null);
        pre.setTool("select");
      }
      if (e.key === "Delete" || e.key === "Backspace") pre.deleteSelected();
      if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === "z") {
        e.preventDefault();
        if (e.shiftKey) pre.redo();
        else pre.undo();
      }
      if (e.key === "Enter" && pre.draft?.tool === "polygon") {
        pre.addPolygon(pre.draft.points);
        pre.setDraft(null);
      }
      if (e.key === "f" || e.key === "F") fit("iso");
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [fit]);

  return (
    <div ref={wrapRef} className="relative h-full min-h-0 w-full bg-viewport">
      <canvas
        ref={canvasRef}
        className="block h-full w-full touch-none"
        onPointerDown={onDown}
        onPointerMove={onMove}
        onPointerUp={onUp}
        onPointerCancel={onUp}
        onContextMenu={(e) => e.preventDefault()}
      />
      <div className="absolute right-3 top-3 flex gap-1">
        <button
          type="button"
          onClick={() => fit("iso")}
          className="h-8 rounded-md border border-border bg-bg/90 px-2.5 text-xs text-muted hover:text-fg"
        >
          Iso
        </button>
        <button
          type="button"
          onClick={() => fit("top")}
          className="h-8 rounded-md border border-border bg-bg/90 px-2.5 text-xs text-muted hover:text-fg"
        >
          Oben
        </button>
        <button
          type="button"
          onClick={() => fit(dim === "2d" ? "top" : "iso")}
          className="h-8 rounded-md border border-border bg-bg/90 px-2.5 text-xs text-muted hover:text-fg"
        >
          Einpassen
        </button>
      </div>
    </div>
  );
}
