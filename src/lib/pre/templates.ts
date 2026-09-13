import { nid } from "./id";
import { defaultMaterials } from "./materials";
import { meshShapes } from "./mesh";
import type { Load, Material, Mesh, Restraint, Shape } from "./types";

export type Snapshot = {
  name: string;
  shapes: Shape[];
  materials: Material[];
  mesh: Mesh | null;
  meshSize: number;
  restraints: Restraint[];
  loads: Load[];
};

export function emptySnapshot(): Snapshot {
  return {
    name: "Modell",
    shapes: [],
    materials: defaultMaterials(),
    mesh: null,
    meshSize: 5,
    restraints: [],
    loads: [],
  };
}

export const TEMPLATE_META = [
  { id: "cantilever", name: "Kragarm" },
  { id: "plate", name: "Gelochte Platte" },
  { id: "bar", name: "Zugstab" },
] as const;

export function templateById(id: string): Snapshot {
  if (id === "plate") return templatePlate();
  if (id === "bar") return templateBar();
  return templateCantilever();
}

function steel(): Material[] {
  return defaultMaterials();
}

export function templateCantilever(): Snapshot {
  const materials = steel();
  const shape: Shape = {
    id: nid(),
    kind: "rect",
    x: 0,
    y: 0,
    w: 100,
    h: 20,
    materialId: "steel",
    holes: [],
  };
  const meshSize = 5;
  const mesh = meshShapes([shape], meshSize);
  const restraints: Restraint[] = [{ id: nid(), target: { type: "edge", shapeId: shape.id, edge: "left" }, ux: true, uy: true }];
  const loads: Load[] = [{ id: nid(), kind: "force", target: { type: "edge", shapeId: shape.id, edge: "right" }, fx: 0, fy: -100 }];
  return { name: "Kragarm", shapes: [shape], materials, mesh, meshSize, restraints, loads };
}

export function templatePlate(): Snapshot {
  const materials = steel();
  const shape: Shape = {
    id: nid(),
    kind: "rect",
    x: 0,
    y: 0,
    w: 200,
    h: 100,
    materialId: "steel",
    holes: [{ id: nid(), cx: 100, cy: 50, r: 20 }],
  };
  const meshSize = 8;
  const mesh = meshShapes([shape], meshSize);
  const restraints: Restraint[] = [
    { id: nid(), target: { type: "edge", shapeId: shape.id, edge: "left" }, ux: true, uy: true },
    { id: nid(), target: { type: "edge", shapeId: shape.id, edge: "right" }, ux: false, uy: true },
  ];
  const loads: Load[] = [{ id: nid(), kind: "force", target: { type: "edge", shapeId: shape.id, edge: "top" }, fx: 0, fy: -200 }];
  return { name: "Gelochte Platte", shapes: [shape], materials, mesh, meshSize, restraints, loads };
}

export function templateBar(): Snapshot {
  const materials = steel();
  const shape: Shape = {
    id: nid(),
    kind: "rect",
    x: 0,
    y: 0,
    w: 80,
    h: 20,
    materialId: "steel",
    holes: [],
  };
  const meshSize = 5;
  const mesh = meshShapes([shape], meshSize);
  const restraints: Restraint[] = [{ id: nid(), target: { type: "edge", shapeId: shape.id, edge: "left" }, ux: true, uy: true }];
  const loads: Load[] = [{ id: nid(), kind: "force", target: { type: "edge", shapeId: shape.id, edge: "right" }, fx: 1000, fy: 0 }];
  return { name: "Zugstab", shapes: [shape], materials, mesh, meshSize, restraints, loads };
}
