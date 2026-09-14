export type Vec2 = { x: number; y: number };
export type Vec3 = { x: number; y: number; z: number };

export type Dim = "2d" | "3d";

export type EdgeName = "left" | "right" | "top" | "bottom" | "boundary";
export type FaceName = "xmin" | "xmax" | "ymin" | "ymax" | "zmin" | "zmax" | "lateral";

export type Hole = { id: string; cx: number; cy: number; r: number };

export type RectShape = {
  id: string;
  kind: "rect";
  x: number;
  y: number;
  w: number;
  h: number;
  depth: number;
  materialId: string;
  holes: Hole[];
};

export type CircleShape = {
  id: string;
  kind: "circle";
  cx: number;
  cy: number;
  r: number;
  depth: number;
  materialId: string;
};

export type PolygonShape = {
  id: string;
  kind: "polygon";
  points: Vec2[];
  depth: number;
  holes: Hole[];
  materialId: string;
};

export type BoxShape = {
  id: string;
  kind: "box";
  x: number;
  y: number;
  z: number;
  w: number;
  h: number;
  d: number;
  materialId: string;
};

export type CylinderShape = {
  id: string;
  kind: "cylinder";
  cx: number;
  cy: number;
  cz: number;
  r: number;
  height: number;
  axis: "x" | "y" | "z";
  materialId: string;
};

export type SphereShape = {
  id: string;
  kind: "sphere";
  cx: number;
  cy: number;
  cz: number;
  r: number;
  materialId: string;
};

export type Shape = RectShape | CircleShape | PolygonShape | BoxShape | CylinderShape | SphereShape;

export type MeshNode = { id: number; x: number; y: number; z: number };

export type ElemType = "CPS3" | "CPS4" | "C3D8" | "C3D6" | "C3D4";

export type MeshElement = {
  id: number;
  type: ElemType;
  nodes: number[];
  materialId: string;
  shapeId: string;
};

export type Mesh = {
  nodes: MeshNode[];
  elements: MeshElement[];
};

export type Material = {
  id: string;
  name: string;
  E: number;
  nu: number;
  density: number;
  thickness: number;
};

export type RestraintTarget =
  | { type: "edge"; shapeId: string; edge: EdgeName }
  | { type: "face"; shapeId: string; face: FaceName }
  | { type: "nodes"; nodeIds: number[] };

export type Restraint = {
  id: string;
  target: RestraintTarget;
  ux: boolean;
  uy: boolean;
  uz: boolean;
};

export type Load =
  | {
      id: string;
      kind: "force";
      target: RestraintTarget;
      fx: number;
      fy: number;
      fz: number;
    }
  | { id: string; kind: "gravity"; g: number };

export type Tool =
  | "select"
  | "rect"
  | "circle"
  | "polygon"
  | "hole"
  | "node"
  | "box"
  | "cylinder"
  | "sphere"
  | "face";

export type Issue = {
  level: "error" | "warn" | "ok";
  code: string;
  message: string;
};

export type MeshQuality = {
  nnode: number;
  nelem: number;
  minAngle: number;
  maxAspect: number;
  nBad: number;
  types: string[];
};

export type SelectedFace = { shapeId: string; face: FaceName };
