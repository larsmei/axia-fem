import { edgeNodes } from "./geometry";
import type { Load, Material, Mesh, Restraint, Shape } from "./types";

function fmt(n: number): string {
  if (!Number.isFinite(n)) return "0";
  const a = Math.abs(n);
  if (a !== 0 && (a >= 1e5 || a < 1e-4)) return n.toExponential(8);
  return (Math.round(n * 1e8) / 1e8).toString();
}

export function resolveNodes(mesh: Mesh, shapes: Shape[], target: Restraint["target"]): number[] {
  if (target.type === "nodes") return [...new Set(target.nodeIds)];
  const shape = shapes.find((s) => s.id === target.shapeId);
  if (!shape) return [];
  const b = shape.kind === "rect" ? Math.min(shape.w, shape.h) : shape.kind === "circle" ? shape.r : 1;
  const tol = Math.max(1e-3, Math.abs(b) * 1e-3);
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
}): string {
  const { name, shapes, mesh, materials, restraints, loads } = opts;
  const used = new Set(mesh.elements.map((e) => e.materialId));
  const mats = materials.filter((m) => used.has(m.id));
  const lines: string[] = [];
  lines.push(`*HEADING`);
  lines.push(`Axia Präprozessor · ${name.replace(/\n/g, " ")}`);
  lines.push(`** units: mm, N, tonne, s → stress in MPa`);
  lines.push(`*NODE`);
  for (const n of mesh.nodes) lines.push(`${n.id}, ${fmt(n.x)}, ${fmt(n.y)}, 0`);

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
    lines.push(`${fmt(m.thickness)}`);
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
    if (r.ux && r.uy) {
      lines.push(`*BOUNDARY`, `${set}, 1, 2`);
    } else if (r.ux) {
      lines.push(`*BOUNDARY`, `${set}, 1, 1`);
    } else if (r.uy) {
      lines.push(`*BOUNDARY`, `${set}, 2, 2`);
    }
  });

  lines.push(`*STEP`);
  lines.push(`*STATIC`);

  loads.forEach((ld, i) => {
    if (ld.kind === "gravity") {
      const mag = Math.abs(ld.g) || 9810;
      lines.push(`*DLOAD`);
      lines.push(`EALL, GRAV, ${fmt(mag)}, 0, -1, 0`);
      return;
    }
    const ids = resolveNodes(mesh, shapes, ld.target);
    if (!ids.length) return;
    const fx = ld.fx / ids.length;
    const fy = ld.fy / ids.length;
    const set = `LOAD_${i + 1}`;
    dumpSet(lines, set, ids);
    lines.push(`*CLOAD`);
    if (Math.abs(ld.fx) > 1e-18) lines.push(`${set}, 1, ${fmt(fx)}`);
    if (Math.abs(ld.fy) > 1e-18) lines.push(`${set}, 2, ${fmt(fy)}`);
  });

  lines.push(`*NODE FILE`);
  lines.push(`U`);
  lines.push(`*EL FILE`);
  lines.push(`S`);
  lines.push(`*END STEP`);
  return lines.join("\n") + "\n";
}
