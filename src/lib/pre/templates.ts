import { nid } from "./id";
import { defaultMaterials } from "./materials";
import { meshShapes } from "./mesh";
import type { Dim, Load, Material, Mesh, Restraint, Shape } from "./types";

export type Snapshot = {
  name: string;
  dim: Dim;
  shapes: Shape[];
  materials: Material[];
  mesh: Mesh | null;
  meshSize: number;
  defaultDepth: number;
  restraints: Restraint[];
  loads: Load[];
};

export function emptySnapshot(): Snapshot {
  return {
    name: "Modell",
    dim: "3d",
    shapes: [],
    materials: defaultMaterials(),
    mesh: null,
    meshSize: 8,
    defaultDepth: 20,
    restraints: [],
    loads: [],
  };
}

export const TEMPLATE_META = [
  { id: "cantilever", name: "Kragarm", dim: "3d" as const },
  { id: "bar", name: "Zugstab", dim: "3d" as const },
  { id: "block", name: "Druckwürfel", dim: "3d" as const },
  { id: "cylinder", name: "Zylinder", dim: "3d" as const },
  { id: "plate", name: "Gelochte Platte", dim: "3d" as const },
  { id: "cantilever2d", name: "Kragarm 2D", dim: "2d" as const },
  { id: "plate2d", name: "Platte 2D", dim: "2d" as const },
  { id: "bar2d", name: "Zugstab 2D", dim: "2d" as const },
] as const;

export function templateById(id: string): Snapshot {
  if (id === "bar") return templateBar3d();
  if (id === "block") return templateBlock();
  if (id === "cylinder") return templateCylinder();
  if (id === "plate") return templatePlate3d();
  if (id === "cantilever2d") return templateCantilever2d();
  if (id === "plate2d") return templatePlate2d();
  if (id === "bar2d") return templateBar2d();
  return templateCantilever3d();
}

function steel(): Material[] {
  return defaultMaterials();
}

function pack(
  name: string,
  dim: Dim,
  shapes: Shape[],
  meshSize: number,
  restraints: Restraint[],
  loads: Load[],
  defaultDepth = 20,
): Snapshot {
  const mesh = meshShapes(shapes, meshSize, dim);
  return { name, dim, shapes, materials: steel(), mesh, meshSize, defaultDepth, restraints, loads };
}

export function templateCantilever(): Snapshot {
  return templateCantilever3d();
}

export function templateCantilever3d(): Snapshot {
  const shape: Shape = {
    id: nid(),
    kind: "box",
    x: 0,
    y: 0,
    z: 0,
    w: 100,
    h: 10,
    d: 20,
    materialId: "steel",
  };
  return pack(
    "Kragarm 3D",
    "3d",
    [shape],
    8,
    [{ id: nid(), target: { type: "face", shapeId: shape.id, face: "xmin" }, ux: true, uy: true, uz: true }],
    [{ id: nid(), kind: "force", target: { type: "face", shapeId: shape.id, face: "xmax" }, fx: 0, fy: 0, fz: -200 }],
    20,
  );
}

export function templateBar3d(): Snapshot {
  const shape: Shape = {
    id: nid(),
    kind: "box",
    x: 0,
    y: 0,
    z: 0,
    w: 80,
    h: 20,
    d: 20,
    materialId: "steel",
  };
  return pack(
    "Zugstab 3D",
    "3d",
    [shape],
    8,
    [{ id: nid(), target: { type: "face", shapeId: shape.id, face: "xmin" }, ux: true, uy: true, uz: true }],
    [{ id: nid(), kind: "force", target: { type: "face", shapeId: shape.id, face: "xmax" }, fx: 2000, fy: 0, fz: 0 }],
    20,
  );
}

export function templateBlock(): Snapshot {
  const shape: Shape = {
    id: nid(),
    kind: "box",
    x: 0,
    y: 0,
    z: 0,
    w: 30,
    h: 30,
    d: 30,
    materialId: "steel",
  };
  return pack(
    "Druckwürfel",
    "3d",
    [shape],
    6,
    [{ id: nid(), target: { type: "face", shapeId: shape.id, face: "zmin" }, ux: true, uy: true, uz: true }],
    [{ id: nid(), kind: "force", target: { type: "face", shapeId: shape.id, face: "zmax" }, fx: 0, fy: 0, fz: -3000 }],
    30,
  );
}

export function templateCylinder(): Snapshot {
  const shape: Shape = {
    id: nid(),
    kind: "cylinder",
    cx: 0,
    cy: 0,
    cz: 30,
    r: 12,
    height: 60,
    axis: "z",
    materialId: "steel",
  };
  return pack(
    "Zylinder",
    "3d",
    [shape],
    6,
    [{ id: nid(), target: { type: "face", shapeId: shape.id, face: "zmin" }, ux: true, uy: true, uz: true }],
    [{ id: nid(), kind: "force", target: { type: "face", shapeId: shape.id, face: "zmax" }, fx: 0, fy: 0, fz: -1500 }],
    60,
  );
}

export function templatePlate3d(): Snapshot {
  const shape: Shape = {
    id: nid(),
    kind: "rect",
    x: 0,
    y: 0,
    w: 120,
    h: 80,
    depth: 8,
    materialId: "steel",
    holes: [{ id: nid(), cx: 60, cy: 40, r: 16 }],
  };
  return pack(
    "Gelochte Platte 3D",
    "3d",
    [shape],
    8,
    [{ id: nid(), target: { type: "face", shapeId: shape.id, face: "xmin" }, ux: true, uy: true, uz: true }],
    [{ id: nid(), kind: "force", target: { type: "face", shapeId: shape.id, face: "xmax" }, fx: 800, fy: 0, fz: 0 }],
    8,
  );
}

export function templateCantilever2d(): Snapshot {
  const shape: Shape = {
    id: nid(),
    kind: "rect",
    x: 0,
    y: 0,
    w: 100,
    h: 20,
    depth: 1,
    materialId: "steel",
    holes: [],
  };
  return pack(
    "Kragarm",
    "2d",
    [shape],
    5,
    [{ id: nid(), target: { type: "edge", shapeId: shape.id, edge: "left" }, ux: true, uy: true, uz: false }],
    [{ id: nid(), kind: "force", target: { type: "edge", shapeId: shape.id, edge: "right" }, fx: 0, fy: -100, fz: 0 }],
    1,
  );
}

export function templatePlate2d(): Snapshot {
  const shape: Shape = {
    id: nid(),
    kind: "rect",
    x: 0,
    y: 0,
    w: 200,
    h: 100,
    depth: 1,
    materialId: "steel",
    holes: [{ id: nid(), cx: 100, cy: 50, r: 20 }],
  };
  return pack(
    "Gelochte Platte",
    "2d",
    [shape],
    8,
    [
      { id: nid(), target: { type: "edge", shapeId: shape.id, edge: "left" }, ux: true, uy: true, uz: false },
      { id: nid(), target: { type: "edge", shapeId: shape.id, edge: "right" }, ux: false, uy: true, uz: false },
    ],
    [{ id: nid(), kind: "force", target: { type: "edge", shapeId: shape.id, edge: "top" }, fx: 0, fy: -200, fz: 0 }],
    1,
  );
}

export function templateBar2d(): Snapshot {
  const shape: Shape = {
    id: nid(),
    kind: "rect",
    x: 0,
    y: 0,
    w: 80,
    h: 20,
    depth: 1,
    materialId: "steel",
    holes: [],
  };
  return pack(
    "Zugstab",
    "2d",
    [shape],
    5,
    [{ id: nid(), target: { type: "edge", shapeId: shape.id, edge: "left" }, ux: true, uy: true, uz: false }],
    [{ id: nid(), kind: "force", target: { type: "edge", shapeId: shape.id, edge: "right" }, fx: 1000, fy: 0, fz: 0 }],
    1,
  );
}

export function templatePlate(): Snapshot {
  return templatePlate3d();
}

export function templateBar(): Snapshot {
  return templateBar3d();
}
