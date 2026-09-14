import { edgeNodes, faceNodes, shapeBounds3 } from "./geometry";
import type { Dim, Load, Material, Mesh, Restraint, Shape } from "./types";

function fmt(n: number): string {
  if (!Number.isFinite(n)) return "0";
  const a = Math.abs(n);
  if (a !== 0 && (a >= 1e5 || a < 1e-4)) return n.toExponential(8);
  return (Math.round(n * 1e8) / 1e8).toString();
}

function isSolidMesh(mesh: Mesh): boolean {
  return mesh.elements.some((e) => e.type.startsWith("C3D"));
}

export function resolveNodes(mesh: Mesh, shapes: Shape[], target: Restraint["target"]): number[] {
  if (target.type === "nodes") return [...new Set(target.nodeIds)];
  const shape = shapes.find((s) => s.id === target.shapeId);
  if (!shape) return [];
  const b = shapeBounds3(shape);
  const span = Math.max(b.w, b.h, b.d, 1);
  const tol = Math.max(1e-3, span * 1e-3);
  if (target.type === "face") return faceNodes(mesh, shape, target.face, tol).map((n) => n.id);
  return edgeNodes(mesh, shape, target.edge, tol).map((n) => n.id);
}

function dumpSet(lines: string[], name: string, ids: number[]) {
  lines.push(`*NSET, NSET=${name}`);
  for (let k = 0; k < ids.length; k += 16) lines.push(ids.slice(k, k + 16).join(", "));
}

export function exportInp(opts: {
  name: string;
  shapes: Shape[];
  mesh: Mesh;
  materials: Material[];
  restraints: Restraint[];
  loads: Load[];
  dim?: Dim;
}): string {
  const { name, shapes, mesh, materials, restraints, loads } = opts;
  const solid = opts.dim === "3d" || isSolidMesh(mesh);
  const used = new Set(mesh.elements.map((e) => e.materialId));
  const mats = materials.filter((m) => used.has(m.id));
  const lines: string[] = [];
  lines.push(`*HEADING`);
  lines.push(`Axia Präprozessor · ${solid ? "3D C3D" : "2D CPS"} · ${name.replace(/\n/g, " ")}`);
  lines.push(`** units: mm, N, tonne, s → stress in MPa`);
  lines.push(`*NODE`);
  for (const n of mesh.nodes) lines.push(`${n.id}, ${fmt(n.x)}, ${fmt(n.y)}, ${fmt(n.z ?? 0)}`);

  const byTypeMat = new Map<string, typeof mesh.elements>();
  for (const el of mesh.elements) {
    const key = `${el.type}|${el.materialId}`;
    const arr = byTypeMat.get(key) ?? [];
    arr.push(el);
    byTypeMat.set(key, arr);
  }
  const elsets = new Set<string>();
  for (const [key, els] of byTypeMat) {
    const [type, matId] = key.split("|");
    const set = `E_${matId.toUpperCase().replace(/[^A-Z0-9_]/g, "_")}`;
    elsets.add(set);
    lines.push(`*ELEMENT, TYPE=${type}, ELSET=${set}`);
    for (const el of els) lines.push(`${el.id}, ${el.nodes.join(", ")}`);
  }

  for (const m of mats) {
    const set = `E_${m.id.toUpperCase().replace(/[^A-Z0-9_]/g, "_")}`;
    if (!elsets.has(set)) continue;
    lines.push(`*SOLID SECTION, ELSET=${set}, MATERIAL=${m.id.toUpperCase()}`);
    if (!solid) lines.push(`${fmt(m.thickness)}`);
    lines.push(`*MATERIAL, NAME=${m.id.toUpperCase()}`);
    lines.push(`*ELASTIC`);
    lines.push(`${fmt(m.E)}, ${fmt(m.nu)}`);
    lines.push(`*DENSITY`);
    lines.push(`${fmt(m.density)}`);
  }

  restraints.forEach((r, i) => {
    const ids = resolveNodes(mesh, shapes, r.target);
    if (!ids.length) return;
    const set = `FIX_${i + 1}`;
    dumpSet(lines, set, ids);
    const ux = r.ux,
      uy = r.uy,
      uz = solid ? r.uz : false;
    if (ux && uy && (uz || !solid)) {
      lines.push(`*BOUNDARY`, `${set}, 1, ${solid ? 3 : 2}`);
    } else {
      lines.push(`*BOUNDARY`);
      if (ux) lines.push(`${set}, 1, 1`);
      if (uy) lines.push(`${set}, 2, 2`);
      if (uz && solid) lines.push(`${set}, 3, 3`);
    }
  });

  lines.push(`*STEP`);
  lines.push(`*STATIC`);

  loads.forEach((ld, i) => {
    if (ld.kind === "gravity") {
      const mag = Math.abs(ld.g) || 9810;
      lines.push(`*DLOAD`);
      if (solid) lines.push(`EALL, GRAV, ${fmt(mag)}, 0, 0, -1`);
      else lines.push(`EALL, GRAV, ${fmt(mag)}, 0, -1, 0`);
      return;
    }
    const ids = resolveNodes(mesh, shapes, ld.target);
    if (!ids.length) return;
    const fx = ld.fx / ids.length;
    const fy = ld.fy / ids.length;
    const fz = (ld.fz ?? 0) / ids.length;
    const set = `LOAD_${i + 1}`;
    dumpSet(lines, set, ids);
    lines.push(`*CLOAD`);
    if (Math.abs(ld.fx) > 1e-18) lines.push(`${set}, 1, ${fmt(fx)}`);
    if (Math.abs(ld.fy) > 1e-18) lines.push(`${set}, 2, ${fmt(fy)}`);
    if (solid && Math.abs(ld.fz) > 1e-18) lines.push(`${set}, 3, ${fmt(fz)}`);
  });

  lines.push(`*NODE FILE`);
  lines.push(`U`);
  lines.push(`*EL FILE`);
  lines.push(`S`);
  lines.push(`*END STEP`);
  return lines.join("\n") + "\n";
}
