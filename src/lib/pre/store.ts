import { create } from "zustand";
import { persist } from "zustand/middleware";
import { checkModel } from "./check";
import { nid } from "./id";
import { meshShapes } from "./mesh";
import { emptySnapshot, templateById, templateCantilever, type Snapshot } from "./templates";
import type { EdgeName, Hole, Issue, Load, Material, Mesh, Restraint, Shape, Tool, Vec2 } from "./types";

type Draft =
  | { tool: "rect"; a: Vec2; b: Vec2 }
  | { tool: "circle"; c: Vec2; r: number }
  | { tool: "polygon"; points: Vec2[] }
  | { tool: "hole"; c: Vec2; r: number }
  | null;

type PreState = Snapshot & {
  selectedShapeId: string | null;
  selectedNodeIds: number[];
  tool: Tool;
  draft: Draft;
  past: Snapshot[];
  future: Snapshot[];
  issues: Issue[];
  setName: (name: string) => void;
  setTool: (t: Tool) => void;
  setDraft: (d: Draft) => void;
  setMeshSize: (n: number) => void;
  selectShape: (id: string | null) => void;
  toggleNode: (id: number, additive?: boolean) => void;
  addRect: (x: number, y: number, w: number, h: number) => void;
  addCircle: (cx: number, cy: number, r: number) => void;
  addPolygon: (points: Vec2[]) => void;
  addHole: (cx: number, cy: number, r: number) => void;
  updateShape: (id: string, patch: Partial<Shape>) => void;
  deleteSelected: () => void;
  setShapeMaterial: (shapeId: string, materialId: string) => void;
  updateMaterial: (id: string, patch: Partial<Material>) => void;
  addCustomMaterial: () => void;
  generateMesh: () => void;
  addEdgeRestraint: (edge: EdgeName, ux: boolean, uy: boolean) => void;
  addNodeRestraint: (ux: boolean, uy: boolean) => void;
  addEdgeForce: (edge: "left" | "right" | "top" | "bottom", fx: number, fy: number) => void;
  addNodeForce: (fx: number, fy: number) => void;
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
    shapes: s.shapes,
    materials: s.materials,
    mesh: s.mesh,
    meshSize: s.meshSize,
    restraints: s.restraints,
    loads: s.loads,
  };
}

const empty = templateCantilever();

export const usePre = create<PreState>()(
  persist(
    (set, get) => {
      const commit = (patch: Partial<Snapshot> & Partial<Pick<PreState, "selectedShapeId" | "selectedNodeIds" | "draft">>) => {
        const cur = get();
        const past = [...cur.past, snapOf(cur)].slice(-40);
        set({ ...patch, past, future: [], draft: patch.draft === undefined ? null : patch.draft });
        queueMicrotask(() => get().refreshIssues());
      };

      return {
        ...empty,
        selectedShapeId: null,
        selectedNodeIds: [],
        tool: "select",
        draft: null,
        past: [],
        future: [],
        issues: [],
        setName: (name) => set({ name }),
        setTool: (tool) => set({ tool, draft: null }),
        setDraft: (draft) => set({ draft }),
        setMeshSize: (meshSize) => set({ meshSize }),
        selectShape: (id) => set({ selectedShapeId: id, selectedNodeIds: [] }),
        toggleNode: (id, additive) =>
          set((s) => {
            if (!additive) return { selectedNodeIds: [id], selectedShapeId: null };
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
          const shape: Shape = {
            id: nid(),
            kind: "rect",
            x,
            y,
            w,
            h,
            materialId: get().materials[0]?.id ?? "steel",
            holes: [],
          };
          commit({ shapes: [...get().shapes, shape], mesh: null, selectedShapeId: shape.id });
        },
        addCircle: (cx, cy, r) => {
          if (r < 1e-6) return;
          const shape: Shape = {
            id: nid(),
            kind: "circle",
            cx,
            cy,
            r,
            materialId: get().materials[0]?.id ?? "steel",
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
            materialId: get().materials[0]?.id ?? "steel",
          };
          commit({ shapes: [...get().shapes, shape], mesh: null, selectedShapeId: shape.id });
        },
        addHole: (cx, cy, r) => {
          const id = get().selectedShapeId;
          const shape = get().shapes.find((s) => s.id === id);
          if (!shape || shape.kind === "circle" || r < 1e-6) return;
          const hole: Hole = { id: nid(), cx, cy, r };
          commit({
            shapes: get().shapes.map((s) => (s.id === id && s.kind !== "circle" ? { ...s, holes: [...s.holes, hole] } : s)),
            mesh: null,
          });
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
              restraints: get().restraints.filter((r) => !(r.target.type === "edge" && r.target.shapeId === selectedShapeId)),
              loads: get().loads.filter((l) => !(l.kind === "force" && l.target.type === "edge" && l.target.shapeId === selectedShapeId)),
              mesh: null,
              selectedShapeId: null,
            });
            return;
          }
          if (selectedNodeIds.length) {
            set({ selectedNodeIds: [] });
          }
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
          const { shapes, meshSize } = get();
          if (!shapes.length) return;
          const mesh = meshShapes(shapes, meshSize);
          commit({ mesh, selectedNodeIds: [] });
        },
        addEdgeRestraint: (edge, ux, uy) => {
          const shapeId = get().selectedShapeId ?? get().shapes[0]?.id;
          if (!shapeId) return;
          const r: Restraint = { id: nid(), target: { type: "edge", shapeId, edge }, ux, uy };
          commit({ restraints: [...get().restraints, r] });
        },
        addNodeRestraint: (ux, uy) => {
          const ids = get().selectedNodeIds;
          if (!ids.length) return;
          commit({
            restraints: [...get().restraints, { id: nid(), target: { type: "nodes", nodeIds: ids }, ux, uy }],
            selectedNodeIds: [],
          });
        },
        addEdgeForce: (edge, fx, fy) => {
          const shapeId = get().selectedShapeId ?? get().shapes[0]?.id;
          if (!shapeId) return;
          commit({
            loads: [...get().loads, { id: nid(), kind: "force", target: { type: "edge", shapeId, edge }, fx, fy }],
          });
        },
        addNodeForce: (fx, fy) => {
          const ids = get().selectedNodeIds;
          if (!ids.length) return;
          commit({
            loads: [...get().loads, { id: nid(), kind: "force", target: { type: "nodes", nodeIds: ids }, fx, fy }],
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
          });
        },
        reset: () => {
          const t = emptySnapshot();
          commit({ ...t, selectedShapeId: null, selectedNodeIds: [] });
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
            }),
          });
        },
      };
    },
    {
      name: "axia-pre-v1",
      partialize: (s) => ({
        name: s.name,
        shapes: s.shapes,
        materials: s.materials,
        mesh: s.mesh,
        meshSize: s.meshSize,
        restraints: s.restraints,
        loads: s.loads,
      }),
    },
  ),
);
