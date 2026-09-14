import { create } from "zustand";
import { persist } from "zustand/middleware";
import { checkModel } from "./check";
import { nid } from "./id";
import { meshShapes } from "./mesh";
import { emptySnapshot, templateById, templateCantilever, type Snapshot } from "./templates";
import type {
  Dim,
  EdgeName,
  FaceName,
  Hole,
  Issue,
  Load,
  Material,
  Mesh,
  Restraint,
  SelectedFace,
  Shape,
  Tool,
  Vec2,
} from "./types";

type Draft =
  | { tool: "rect"; a: Vec2; b: Vec2 }
  | { tool: "circle"; c: Vec2; r: number }
  | { tool: "polygon"; points: Vec2[] }
  | { tool: "hole"; c: Vec2; r: number }
  | { tool: "box"; a: Vec2; b: Vec2 }
  | { tool: "cylinder"; c: Vec2; r: number }
  | { tool: "sphere"; c: Vec2; r: number }
  | null;

type PreState = Snapshot & {
  selectedShapeId: string | null;
  selectedNodeIds: number[];
  selectedFace: SelectedFace | null;
  tool: Tool;
  draft: Draft;
  past: Snapshot[];
  future: Snapshot[];
  issues: Issue[];
  setName: (name: string) => void;
  setDim: (dim: Dim) => void;
  setTool: (t: Tool) => void;
  setDraft: (d: Draft) => void;
  setMeshSize: (n: number) => void;
  setDefaultDepth: (n: number) => void;
  selectShape: (id: string | null) => void;
  selectFace: (f: SelectedFace | null) => void;
  toggleNode: (id: number, additive?: boolean) => void;
  addRect: (x: number, y: number, w: number, h: number) => void;
  addCircle: (cx: number, cy: number, r: number) => void;
  addPolygon: (points: Vec2[]) => void;
  addHole: (cx: number, cy: number, r: number) => void;
  addBox: (x: number, y: number, z: number, w: number, h: number, d: number) => void;
  addCylinder: (cx: number, cy: number, cz: number, r: number, height: number, axis?: "x" | "y" | "z") => void;
  addSphere: (cx: number, cy: number, cz: number, r: number) => void;
  updateShape: (id: string, patch: Partial<Shape>) => void;
  deleteSelected: () => void;
  setShapeMaterial: (shapeId: string, materialId: string) => void;
  updateMaterial: (id: string, patch: Partial<Material>) => void;
  addCustomMaterial: () => void;
  generateMesh: () => void;
  addEdgeRestraint: (edge: EdgeName, ux: boolean, uy: boolean, uz?: boolean) => void;
  addFaceRestraint: (face: FaceName, ux: boolean, uy: boolean, uz: boolean) => void;
  addNodeRestraint: (ux: boolean, uy: boolean, uz?: boolean) => void;
  addEdgeForce: (edge: "left" | "right" | "top" | "bottom", fx: number, fy: number, fz?: number) => void;
  addFaceForce: (face: FaceName, fx: number, fy: number, fz: number) => void;
  addNodeForce: (fx: number, fy: number, fz?: number) => void;
  toggleGravity: (on: boolean) => void;
  removeRestraint: (id: string) => void;
  removeLoad: (id: string) => void;
  loadTemplate: (id: string) => void;
  reset: () => void;
  undo: () => void;
  redo: () => void;
  refreshIssues: () => void;
};

function snapOf(s: Snapshot): Snapshot {
  return {
    name: s.name,
    dim: s.dim,
    shapes: s.shapes,
    materials: s.materials,
    mesh: s.mesh,
    meshSize: s.meshSize,
    defaultDepth: s.defaultDepth,
    restraints: s.restraints,
    loads: s.loads,
  };
}

function mat0(get: () => PreState): string {
  return get().materials[0]?.id ?? "steel";
}

const empty = templateCantilever();

export const usePre = create<PreState>()(
  persist(
    (set, get) => {
      const commit = (
        patch: Partial<Snapshot> & Partial<Pick<PreState, "selectedShapeId" | "selectedNodeIds" | "selectedFace" | "draft">>,
      ) => {
        const cur = get();
        const past = [...cur.past, snapOf(cur)].slice(-40);
        set({
          ...patch,
          past,
          future: [],
          draft: patch.draft === undefined ? null : patch.draft,
        });
        queueMicrotask(() => get().refreshIssues());
      };

      return {
        ...empty,
        selectedShapeId: empty.shapes[0]?.id ?? null,
        selectedNodeIds: [],
        selectedFace: null,
        tool: "select",
        draft: null,
        past: [],
        future: [],
        issues: [],
        setName: (name) => set({ name }),
        setDim: (dim) => {
          const cur = get();
          if (cur.dim === dim) return;
          let shapes = cur.shapes;
          if (dim === "2d") {
            shapes = shapes.filter((s) => s.kind === "rect" || s.kind === "circle" || s.kind === "polygon");
          }
          commit({
            dim,
            shapes,
            mesh: null,
            meshSize: dim === "3d" ? Math.max(cur.meshSize, 6) : cur.meshSize,
            selectedFace: null,
          });
        },
        setTool: (tool) => set({ tool, draft: null }),
        setDraft: (draft) => set({ draft }),
        setMeshSize: (meshSize) => set({ meshSize }),
        setDefaultDepth: (defaultDepth) => set({ defaultDepth: Math.max(0.5, defaultDepth) }),
        selectShape: (id) => set({ selectedShapeId: id, selectedNodeIds: [], selectedFace: null }),
        selectFace: (f) => set({ selectedFace: f, selectedShapeId: f?.shapeId ?? get().selectedShapeId, selectedNodeIds: [] }),
        toggleNode: (id, additive) =>
          set((s) => {
            if (!additive) return { selectedNodeIds: [id], selectedShapeId: null, selectedFace: null };
            const has = s.selectedNodeIds.includes(id);
            return { selectedNodeIds: has ? s.selectedNodeIds.filter((n) => n !== id) : [...s.selectedNodeIds, id] };
          }),
        addRect: (x, y, w, h) => {
          if (w < 0) {
            x += w;
            w = -w;
          }
          if (h < 0) {
            y += h;
            h = -h;
          }
          if (w < 1e-6 || h < 1e-6) return;
          const dim = get().dim;
          if (dim === "3d") {
            get().addBox(x, y, 0, w, h, get().defaultDepth);
            return;
          }
          const shape: Shape = {
            id: nid(),
            kind: "rect",
            x,
            y,
            w,
            h,
            depth: get().defaultDepth,
            materialId: mat0(get),
            holes: [],
          };
          commit({ shapes: [...get().shapes, shape], mesh: null, selectedShapeId: shape.id });
        },
        addCircle: (cx, cy, r) => {
          if (r < 1e-6) return;
          const dim = get().dim;
          if (dim === "3d") {
            get().addCylinder(cx, cy, get().defaultDepth / 2, r, get().defaultDepth, "z");
            return;
          }
          const shape: Shape = {
            id: nid(),
            kind: "circle",
            cx,
            cy,
            r,
            depth: get().defaultDepth,
            materialId: mat0(get),
          };
          commit({ shapes: [...get().shapes, shape], mesh: null, selectedShapeId: shape.id });
        },
        addPolygon: (points) => {
          if (points.length < 3) return;
          const shape: Shape = {
            id: nid(),
            kind: "polygon",
            points,
            holes: [],
            depth: get().defaultDepth,
            materialId: mat0(get),
          };
          commit({ shapes: [...get().shapes, shape], mesh: null, selectedShapeId: shape.id });
        },
        addHole: (cx, cy, r) => {
          const id = get().selectedShapeId;
          const shape = get().shapes.find((s) => s.id === id);
          if (!shape || (shape.kind !== "rect" && shape.kind !== "polygon") || r < 1e-6) return;
          const hole: Hole = { id: nid(), cx, cy, r };
          commit({
            shapes: get().shapes.map((s) =>
              s.id === id && (s.kind === "rect" || s.kind === "polygon") ? { ...s, holes: [...s.holes, hole] } : s,
            ),
            mesh: null,
          });
        },
        addBox: (x, y, z, w, h, d) => {
          if (w < 0) {
            x += w;
            w = -w;
          }
          if (h < 0) {
            y += h;
            h = -h;
          }
          if (d < 0) {
            z += d;
            d = -d;
          }
          if (w < 1e-6 || h < 1e-6 || d < 1e-6) return;
          const shape: Shape = { id: nid(), kind: "box", x, y, z, w, h, d, materialId: mat0(get) };
          commit({ shapes: [...get().shapes, shape], mesh: null, selectedShapeId: shape.id, dim: "3d" });
        },
        addCylinder: (cx, cy, cz, r, height, axis = "z") => {
          if (r < 1e-6 || height < 1e-6) return;
          const shape: Shape = { id: nid(), kind: "cylinder", cx, cy, cz, r, height, axis, materialId: mat0(get) };
          commit({ shapes: [...get().shapes, shape], mesh: null, selectedShapeId: shape.id, dim: "3d" });
        },
        addSphere: (cx, cy, cz, r) => {
          if (r < 1e-6) return;
          const shape: Shape = { id: nid(), kind: "sphere", cx, cy, cz, r, materialId: mat0(get) };
          commit({ shapes: [...get().shapes, shape], mesh: null, selectedShapeId: shape.id, dim: "3d" });
        },
        updateShape: (id, patch) => {
          commit({
            shapes: get().shapes.map((s) => (s.id === id ? ({ ...s, ...patch } as Shape) : s)),
            mesh: null,
          });
        },
        deleteSelected: () => {
          const { selectedShapeId, selectedNodeIds } = get();
          if (selectedShapeId) {
            commit({
              shapes: get().shapes.filter((s) => s.id !== selectedShapeId),
              restraints: get().restraints.filter(
                (r) =>
                  !((r.target.type === "edge" || r.target.type === "face") && r.target.shapeId === selectedShapeId),
              ),
              loads: get().loads.filter(
                (l) =>
                  !(
                    l.kind === "force" &&
                    (l.target.type === "edge" || l.target.type === "face") &&
                    l.target.shapeId === selectedShapeId
                  ),
              ),
              mesh: null,
              selectedShapeId: null,
              selectedFace: null,
            });
            return;
          }
          if (selectedNodeIds.length) set({ selectedNodeIds: [] });
        },
        setShapeMaterial: (shapeId, materialId) => {
          commit({
            shapes: get().shapes.map((s) => (s.id === shapeId ? { ...s, materialId } : s)),
            mesh: get().mesh
              ? {
                  ...get().mesh!,
                  elements: get().mesh!.elements.map((e) => (e.shapeId === shapeId ? { ...e, materialId } : e)),
                }
              : null,
          });
        },
        updateMaterial: (id, patch) => {
          commit({ materials: get().materials.map((m) => (m.id === id ? { ...m, ...patch } : m)) });
        },
        addCustomMaterial: () => {
          const m: Material = {
            id: `mat_${nid()}`,
            name: "Custom",
            E: 210000,
            nu: 0.3,
            density: 7.85e-9,
            thickness: 1,
          };
          commit({ materials: [...get().materials, m] });
        },
        generateMesh: () => {
          const { shapes, meshSize, dim } = get();
          if (!shapes.length) return;
          const mesh = meshShapes(shapes, meshSize, dim);
          commit({ mesh, selectedNodeIds: [] });
        },
        addEdgeRestraint: (edge, ux, uy, uz) => {
          const shapeId = get().selectedShapeId ?? get().shapes[0]?.id;
          if (!shapeId) return;
          const solid = get().dim === "3d";
          const r: Restraint = {
            id: nid(),
            target: { type: "edge", shapeId, edge },
            ux,
            uy,
            uz: uz ?? solid,
          };
          commit({ restraints: [...get().restraints, r] });
        },
        addFaceRestraint: (face, ux, uy, uz) => {
          const shapeId = get().selectedFace?.shapeId ?? get().selectedShapeId ?? get().shapes[0]?.id;
          if (!shapeId) return;
          commit({
            restraints: [
              ...get().restraints,
              { id: nid(), target: { type: "face", shapeId, face }, ux, uy, uz },
            ],
          });
        },
        addNodeRestraint: (ux, uy, uz) => {
          const ids = get().selectedNodeIds;
          if (!ids.length) return;
          const solid = get().dim === "3d";
          commit({
            restraints: [
              ...get().restraints,
              { id: nid(), target: { type: "nodes", nodeIds: ids }, ux, uy, uz: uz ?? solid },
            ],
            selectedNodeIds: [],
          });
        },
        addEdgeForce: (edge, fx, fy, fz = 0) => {
          const shapeId = get().selectedShapeId ?? get().shapes[0]?.id;
          if (!shapeId) return;
          commit({
            loads: [...get().loads, { id: nid(), kind: "force", target: { type: "edge", shapeId, edge }, fx, fy, fz }],
          });
        },
        addFaceForce: (face, fx, fy, fz) => {
          const shapeId = get().selectedFace?.shapeId ?? get().selectedShapeId ?? get().shapes[0]?.id;
          if (!shapeId) return;
          commit({
            loads: [...get().loads, { id: nid(), kind: "force", target: { type: "face", shapeId, face }, fx, fy, fz }],
          });
        },
        addNodeForce: (fx, fy, fz = 0) => {
          const ids = get().selectedNodeIds;
          if (!ids.length) return;
          commit({
            loads: [...get().loads, { id: nid(), kind: "force", target: { type: "nodes", nodeIds: ids }, fx, fy, fz }],
            selectedNodeIds: [],
          });
        },
        toggleGravity: (on) => {
          const rest = get().loads.filter((l) => l.kind !== "gravity");
          commit({ loads: on ? [...rest, { id: nid(), kind: "gravity", g: 9810 }] : rest });
        },
        removeRestraint: (id) => commit({ restraints: get().restraints.filter((r) => r.id !== id) }),
        removeLoad: (id) => commit({ loads: get().loads.filter((l) => l.id !== id) }),
        loadTemplate: (id) => {
          const t = templateById(id);
          commit({
            ...t,
            selectedShapeId: t.shapes[0]?.id ?? null,
            selectedNodeIds: [],
            selectedFace: null,
          });
        },
        reset: () => {
          const t = emptySnapshot();
          commit({ ...t, selectedShapeId: null, selectedNodeIds: [], selectedFace: null });
        },
        undo: () => {
          const { past, future } = get();
          if (!past.length) return;
          const prev = past[past.length - 1];
          set({
            ...prev,
            past: past.slice(0, -1),
            future: [snapOf(get()), ...future].slice(0, 40),
            selectedShapeId: null,
            selectedNodeIds: [],
            selectedFace: null,
            draft: null,
          });
          queueMicrotask(() => get().refreshIssues());
        },
        redo: () => {
          const { future } = get();
          if (!future.length) return;
          const next = future[0];
          set({
            ...next,
            past: [...get().past, snapOf(get())].slice(-40),
            future: future.slice(1),
            selectedShapeId: null,
            selectedNodeIds: [],
            selectedFace: null,
            draft: null,
          });
          queueMicrotask(() => get().refreshIssues());
        },
        refreshIssues: () => {
          const s = get();
          set({
            issues: checkModel({
              shapes: s.shapes,
              mesh: s.mesh,
              materials: s.materials,
              restraints: s.restraints,
              loads: s.loads,
              dim: s.dim,
            }),
          });
        },
      };
    },
    {
      name: "axia-pre-v2",
      partialize: (s) => ({
        name: s.name,
        dim: s.dim,
        shapes: s.shapes,
        materials: s.materials,
        mesh: s.mesh,
        meshSize: s.meshSize,
        defaultDepth: s.defaultDepth,
        restraints: s.restraints,
        loads: s.loads,
      }),
    },
  ),
);

export type { Draft, Mesh };
