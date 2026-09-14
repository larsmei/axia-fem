import { resolveNodes } from "./export-inp";
import { meshQuality } from "./mesh";
import type { Dim, Issue, Load, Material, Mesh, Restraint, Shape } from "./types";

export function checkModel(opts: {
  shapes: Shape[];
  mesh: Mesh | null;
  materials: Material[];
  restraints: Restraint[];
  loads: Load[];
  dim?: Dim;
}): Issue[] {
  const { shapes, mesh, materials, restraints, loads } = opts;
  const dim = opts.dim ?? "2d";
  const solid = dim === "3d" || !!mesh?.elements.some((e) => e.type.startsWith("C3D"));
  const issues: Issue[] = [];
  if (!shapes.length) issues.push({ level: "error", code: "geo", message: "Keine Geometrie." });
  if (!mesh || !mesh.elements.length) {
    issues.push({ level: "error", code: "mesh", message: "Kein Netz. Zuerst vernetzen." });
  } else {
    const q = meshQuality(mesh);
    if (q.nBad) {
      issues.push({
        level: "warn",
        code: "quality",
        message: `${q.nBad} Elemente mit spitzem Winkel (<12°) oder Aspekt > 12. min∠ ${q.minAngle.toFixed(1)}°.`,
      });
    }
    if (mesh.nodes.length > 2500) {
      issues.push({
        level: "warn",
        code: "size",
        message: `${mesh.nodes.length} Knoten — Browser-Solver mag kleine Netze.`,
      });
    }
  }
  const usedMats = new Set(mesh?.elements.map((e) => e.materialId) ?? []);
  for (const id of usedMats) {
    const m = materials.find((x) => x.id === id);
    if (!m) issues.push({ level: "error", code: "mat", message: `Material ${id} fehlt.` });
    else {
      if (!(m.E > 0)) issues.push({ level: "error", code: "mat", message: `${m.name}: E muss > 0 sein.` });
      if (!(m.nu > 0 && m.nu < 0.5)) issues.push({ level: "error", code: "mat", message: `${m.name}: ν zwischen 0 und 0,5.` });
      if (!solid && !(m.thickness > 0)) issues.push({ level: "error", code: "mat", message: `${m.name}: Dicke muss > 0 sein.` });
    }
  }
  if (mesh) {
    const fixed = new Set<number>();
    let hasUx = false,
      hasUy = false,
      hasUz = false;
    for (const r of restraints) {
      const ids = resolveNodes(mesh, shapes, r.target);
      if (!ids.length) issues.push({ level: "error", code: "bc", message: "Lager ohne Knoten — Fläche nach Vernetzung prüfen." });
      if (!r.ux && !r.uy && !(solid && r.uz)) issues.push({ level: "warn", code: "bc", message: "Lager ohne gesperrten DOF." });
      ids.forEach((id) => fixed.add(id));
      if (r.ux) hasUx = true;
      if (r.uy) hasUy = true;
      if (r.uz) hasUz = true;
    }
    if (!restraints.length) issues.push({ level: "error", code: "bc", message: "Keine Lager — Mechanismus." });
    else if (solid && !(hasUx && hasUy && hasUz)) {
      issues.push({ level: "warn", code: "bc", message: "3D: Ux, Uy und Uz sollten irgendwo gesperrt sein." });
    }
    const hasLoad = loads.some((l) => l.kind === "gravity" || (l.kind === "force" && (l.fx || l.fy || l.fz)));
    if (!hasLoad) issues.push({ level: "warn", code: "load", message: "Keine Last." });
    for (const l of loads) {
      if (l.kind === "force") {
        const ids = resolveNodes(mesh, shapes, l.target);
        if (!ids.length) issues.push({ level: "error", code: "load", message: "Last ohne Knoten." });
      }
    }
    if (hasLoad && !fixed.size) issues.push({ level: "error", code: "bc", message: "Last ohne Lager." });
  }
  if (!issues.length) issues.push({ level: "ok", code: "ok", message: "Modell vollständig — bereit zum Export." });
  return issues;
}

export function hasErrors(issues: Issue[]): boolean {
  return issues.some((i) => i.level === "error");
}
