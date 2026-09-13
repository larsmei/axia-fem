export type Vec2 = { x: number; y: number };

export type EdgeName = "left" | "right" | "top" | "bottom" | "boundary";

export type Hole = { id: string; cx: number; cy: number; r: number };

export type RectShape = {
  id: string;
  kind: "rect";
  x: number;
  y: number;
  w: number;
  h: number;
  materialId: string;
  holes: Hole[];
};

export type CircleShape = {
  id: string;
  kind: "circle";
  cx: number;
  cy: number;
  r: number;
  materialId: string;
};

export type PolygonShape = {
  id: string;
  kind: "polygon";
  points: Vec2[];
  holes: Hole[];
  materialId: string;
};

export type Shape = RectShape | CircleShape | PolygonShape;

export type MeshNode = { id: number; x: number; y: number };

export type MeshElement = {
  id: number;
  type: "CPS3" | "CPS4";
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
  | { type: "nodes"; nodeIds: number[] };

export type Restraint = {
  id: string;
  target: RestraintTarget;
  ux: boolean;
  uy: boolean;
};

export type Load =
  | {
      id: string;
      kind: "force";
      target: RestraintTarget;
      fx: number;
      fy: number;
    }
  | { id: string; kind: "gravity"; g: number };

export type Tool = "select" | "rect" | "circle" | "polygon" | "hole" | "node";

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
};
