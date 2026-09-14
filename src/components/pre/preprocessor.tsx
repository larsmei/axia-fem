import { useEffect, useMemo, useState, type ReactNode } from "react";
import { useNavigate } from "@tanstack/react-router";
import {
  Box,
  BoxSelect,
  CheckCircle2,
  Circle,
  Cylinder,
  Download,
  Grip,
  Hexagon,
  Layers,
  MousePointer2,
  Play,
  Plus,
  Redo2,
  Square,
  Trash2,
  Undo2,
} from "lucide-react";
import { AppMark, AppNav } from "@/components/app-nav";
import { Button } from "@/components/ui/button";
import { PreViewport } from "@/components/pre/viewport";
import { checkModel, hasErrors } from "@/lib/pre/check";
import { sendInpToSolver } from "@/lib/pre/bridge";
import { exportInp } from "@/lib/pre/export-inp";
import { faceLabel, shapeKindLabel } from "@/lib/pre/geometry";
import { elemTypeSummary, meshQuality } from "@/lib/pre/mesh";
import { usePre } from "@/lib/pre/store";
import { TEMPLATE_META } from "@/lib/pre/templates";
import type { FaceName, Tool } from "@/lib/pre/types";
import { cn, downloadText } from "@/lib/utils";

type Step = "geo" | "mesh" | "mat" | "bc" | "check";

const STEPS: { id: Step; n: string; label: string }[] = [
  { id: "geo", n: "1", label: "Geometrie" },
  { id: "mesh", n: "2", label: "Netz" },
  { id: "mat", n: "3", label: "Material" },
  { id: "bc", n: "4", label: "Lager / Last" },
  { id: "check", n: "5", label: "Prüfen" },
];

const FACES: { id: FaceName; label: string }[] = [
  { id: "xmin", label: "X−" },
  { id: "xmax", label: "X+" },
  { id: "ymin", label: "Y−" },
  { id: "ymax", label: "Y+" },
  { id: "zmin", label: "Z−" },
  { id: "zmax", label: "Z+" },
];

export function Preprocessor() {
  const navigate = useNavigate();
  const [step, setStep] = useState<Step>("geo");
  const name = usePre((s) => s.name);
  const dim = usePre((s) => s.dim);
  const shapes = usePre((s) => s.shapes);
  const mesh = usePre((s) => s.mesh);
  const issues = usePre((s) => s.issues);
  const refreshIssues = usePre((s) => s.refreshIssues);

  useEffect(() => {
    refreshIssues();
  }, [refreshIssues]);

  const q = useMemo(() => (mesh ? meshQuality(mesh) : null), [mesh]);
  const ready = !hasErrors(issues) && !!mesh;

  const send = () => {
    const s = usePre.getState();
    const list = checkModel(s);
    if (hasErrors(list) || !s.mesh) return;
    const inp = exportInp({
      name: s.name,
      shapes: s.shapes,
      mesh: s.mesh,
      materials: s.materials,
      restraints: s.restraints,
      loads: s.loads,
      dim: s.dim,
    });
    sendInpToSolver(inp, s.name);
    void navigate({ to: "/" });
  };

  const download = () => {
    const s = usePre.getState();
    if (!s.mesh) return;
    const inp = exportInp({
      name: s.name,
      shapes: s.shapes,
      mesh: s.mesh,
      materials: s.materials,
      restraints: s.restraints,
      loads: s.loads,
      dim: s.dim,
    });
    downloadText(`${s.name.replace(/\s+/g, "_").toLowerCase() || "model"}.inp`, inp);
  };

  return (
    <main className="flex h-dvh flex-col overflow-hidden bg-bg text-fg">
      <header className="flex shrink-0 items-center gap-2 border-b border-border bg-surface px-3 py-2 lg:px-4">
        <AppMark />
        <div className="min-w-0 flex-1">
          <div className="flex items-baseline gap-2">
            <h1 className="text-sm font-medium tracking-tight">Axia</h1>
            <p className="hidden text-xs text-muted sm:block">{dim === "3d" ? "Präprozessor · 3D C3D8" : "Präprozessor · 2D CPS"}</p>
          </div>
        </div>
        <AppNav active="pre" />
        <div className="hidden rounded-md border border-border bg-surface-2 p-0.5 sm:flex">
          <button
            type="button"
            onClick={() => usePre.getState().setDim("2d")}
            className={cn("h-7 rounded px-2 text-xs font-medium", dim === "2d" ? "bg-surface text-fg" : "text-muted hover:text-fg")}
          >
            2D
          </button>
          <button
            type="button"
            onClick={() => usePre.getState().setDim("3d")}
            className={cn("h-7 rounded px-2 text-xs font-medium", dim === "3d" ? "bg-surface text-fg" : "text-muted hover:text-fg")}
          >
            3D
          </button>
        </div>
        <input
          value={name}
          onChange={(e) => usePre.getState().setName(e.target.value)}
          className="hidden h-8 w-32 rounded-md border border-border bg-surface-2 px-2 text-xs text-fg lg:block"
          aria-label="Modellname"
        />
        <Button variant="ghost" size="icon-sm" onClick={() => usePre.getState().undo()} aria-label="Rückgängig">
          <Undo2 />
        </Button>
        <Button variant="ghost" size="icon-sm" onClick={() => usePre.getState().redo()} aria-label="Wiederholen">
          <Redo2 />
        </Button>
        <Button variant="secondary" size="sm" onClick={download} disabled={!mesh}>
          <Download />
          <span className="hidden sm:inline">INP</span>
        </Button>
        <Button onClick={send} disabled={!ready} className="min-w-24">
          <Play />
          Lösen
        </Button>
      </header>

      <div className="flex gap-1 overflow-x-auto border-b border-border px-2 py-1">
        {STEPS.map((s) => (
          <button
            key={s.id}
            type="button"
            onClick={() => setStep(s.id)}
            className={cn(
              "flex h-9 shrink-0 items-center gap-2 rounded-md px-2.5 text-xs font-medium",
              step === s.id ? "bg-surface-2 text-fg" : "text-muted hover:text-fg",
            )}
          >
            <span className="font-mono text-subtle">{s.n}</span>
            {s.label}
          </button>
        ))}
      </div>

      <div className="flex min-h-0 flex-1">
        <aside className="hidden w-12 shrink-0 flex-col items-center gap-1 border-r border-border bg-surface py-2 sm:flex">
          <ToolBtn id="select" icon={MousePointer2} label="Auswählen" />
          {dim === "3d" ? (
            <>
              <ToolBtn id="box" icon={Box} label="Quader" />
              <ToolBtn id="cylinder" icon={Cylinder} label="Zylinder" />
              <ToolBtn id="sphere" icon={Circle} label="Kugel" />
              <ToolBtn id="polygon" icon={Hexagon} label="Polygon (extrudiert)" />
            </>
          ) : (
            <>
              <ToolBtn id="rect" icon={Square} label="Rechteck" />
              <ToolBtn id="circle" icon={Circle} label="Kreis" />
              <ToolBtn id="polygon" icon={Hexagon} label="Polygon" />
            </>
          )}
          <ToolBtn id="hole" icon={BoxSelect} label="Loch" />
          <ToolBtn id="node" icon={Grip} label="Knoten" />
          <div className="mt-auto flex flex-col gap-1">
            <Button variant="ghost" size="icon-sm" onClick={() => usePre.getState().deleteSelected()} aria-label="Löschen">
              <Trash2 />
            </Button>
          </div>
        </aside>

        <div className="flex min-h-0 min-w-0 flex-1 flex-col">
          <div className="flex gap-1 overflow-x-auto border-b border-border px-2 py-1 sm:hidden">
            <ToolBtn id="select" icon={MousePointer2} label="Auswahl" wide />
            {dim === "3d" ? (
              <>
                <ToolBtn id="box" icon={Box} label="Quader" wide />
                <ToolBtn id="cylinder" icon={Cylinder} label="Zylinder" wide />
                <ToolBtn id="sphere" icon={Circle} label="Kugel" wide />
              </>
            ) : (
              <>
                <ToolBtn id="rect" icon={Square} label="Rechteck" wide />
                <ToolBtn id="circle" icon={Circle} label="Kreis" wide />
              </>
            )}
            <ToolBtn id="polygon" icon={Hexagon} label="Polygon" wide />
            <ToolBtn id="hole" icon={BoxSelect} label="Loch" wide />
            <ToolBtn id="node" icon={Grip} label="Knoten" wide />
          </div>
          <div className="min-h-0 flex-1">
            <PreViewport />
          </div>
        </div>

        <aside className="hidden w-72 shrink-0 flex-col border-l border-border bg-surface lg:flex">
          <div className="min-h-0 flex-1 overflow-y-auto p-3">
            {step === "geo" && <GeoPanel />}
            {step === "mesh" && <MeshPanel />}
            {step === "mat" && <MatPanel />}
            {step === "bc" && <BcPanel />}
            {step === "check" && <CheckPanel />}
          </div>
        </aside>
      </div>

      <div className="max-h-[34vh] overflow-y-auto border-t border-border bg-surface p-3 lg:hidden">
        {step === "geo" && <GeoPanel />}
        {step === "mesh" && <MeshPanel />}
        {step === "mat" && <MatPanel />}
        {step === "bc" && <BcPanel />}
        {step === "check" && <CheckPanel />}
      </div>

      <footer className="flex shrink-0 items-center justify-between gap-3 border-t border-border bg-surface px-3 py-2 text-xs text-muted">
        <p className="min-w-0 truncate">
          {dim.toUpperCase()} · {shapes.length} Körper · {q ? `${q.nnode} Knoten · ${q.nelem} ${elemTypeSummary(mesh)}` : "kein Netz"} · mm / N / MPa
        </p>
        <p className="hidden shrink-0 sm:block">{issues.find((i) => i.level === "error")?.message ?? "Präprozessor"}</p>
      </footer>
    </main>
  );
}

function ToolBtn({ id, icon: Icon, label, wide }: { id: Tool; icon: typeof Square; label: string; wide?: boolean }) {
  const tool = usePre((s) => s.tool);
  return (
    <button
      type="button"
      title={label}
      onClick={() => usePre.getState().setTool(id)}
      className={cn(
        "inline-flex items-center justify-center rounded-md",
        wide ? "h-9 gap-1.5 px-2.5 text-xs" : "size-9",
        tool === id ? "bg-surface-2 text-fg" : "text-muted hover:text-fg",
      )}
    >
      <Icon className="size-4" />
      {wide && label}
    </button>
  );
}

function GeoPanel() {
  const dim = usePre((s) => s.dim);
  const shapes = usePre((s) => s.shapes);
  const selected = usePre((s) => s.shapes.find((x) => x.id === s.selectedShapeId));
  const defaultDepth = usePre((s) => s.defaultDepth);
  const templates = TEMPLATE_META.filter((t) => t.dim === dim);
  return (
    <div className="flex flex-col gap-3">
      <div className="lg:hidden">
        <Section title="Modus">
          <div className="flex rounded-md border border-border bg-surface-2 p-0.5">
            <button
              type="button"
              onClick={() => usePre.getState().setDim("2d")}
              className={cn("h-8 flex-1 rounded text-xs font-medium", dim === "2d" ? "bg-surface text-fg" : "text-muted")}
            >
              2D Scheibe
            </button>
            <button
              type="button"
              onClick={() => usePre.getState().setDim("3d")}
              className={cn("h-8 flex-1 rounded text-xs font-medium", dim === "3d" ? "bg-surface text-fg" : "text-muted")}
            >
              3D Volumen
            </button>
          </div>
        </Section>
      </div>
      <Section title="Vorlagen">
        <div className="flex flex-wrap gap-1.5">
          {templates.map((t) => (
            <Button key={t.id} variant="secondary" size="sm" onClick={() => usePre.getState().loadTemplate(t.id)}>
              {t.name}
            </Button>
          ))}
          <Button variant="ghost" size="sm" onClick={() => usePre.getState().reset()}>
            Neu
          </Button>
        </div>
      </Section>
      {dim === "3d" && (
        <Section title="Einfügen">
          <div className="flex flex-wrap gap-1.5">
            <Button variant="secondary" size="sm" onClick={() => usePre.getState().addBox(0, 0, 0, 40, 20, 20)}>
              Quader
            </Button>
            <Button variant="secondary" size="sm" onClick={() => usePre.getState().addCylinder(0, 0, 20, 12, 40, "z")}>
              Zylinder
            </Button>
            <Button variant="secondary" size="sm" onClick={() => usePre.getState().addSphere(0, 0, 15, 15)}>
              Kugel
            </Button>
          </div>
          <Num label="Höhe Z" value={defaultDepth} unit="mm" onChange={(v) => usePre.getState().setDefaultDepth(v)} />
        </Section>
      )}
      <Section title="Körper">
        {shapes.length === 0 && (
          <p className="text-xs text-muted">{dim === "3d" ? "Quader, Zylinder oder Kugel zeichnen." : "Rechteck, Kreis oder Polygon zeichnen."}</p>
        )}
        <ul className="flex flex-col gap-1">
          {shapes.map((s, i) => (
            <li key={s.id}>
              <button
                type="button"
                onClick={() => usePre.getState().selectShape(s.id)}
                className={cn(
                  "flex h-8 w-full items-center rounded-md px-2 text-left text-xs",
                  selected?.id === s.id ? "bg-surface-2 text-fg" : "text-muted hover:text-fg",
                )}
              >
                {shapeKindLabel(s)} {i + 1}
              </button>
            </li>
          ))}
        </ul>
      </Section>
      {selected?.kind === "rect" && (
        <Section title="Abmessungen">
          <Num label="x" value={selected.x} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { x: v })} />
          <Num label="y" value={selected.y} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { y: v })} />
          <Num label="b" value={selected.w} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { w: Math.max(0.1, v) })} />
          <Num label="h" value={selected.h} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { h: Math.max(0.1, v) })} />
          {dim === "3d" && (
            <Num label="t" value={selected.depth} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { depth: Math.max(0.1, v) })} />
          )}
        </Section>
      )}
      {selected?.kind === "box" && (
        <Section title="Abmessungen">
          <Num label="x" value={selected.x} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { x: v })} />
          <Num label="y" value={selected.y} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { y: v })} />
          <Num label="z" value={selected.z} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { z: v })} />
          <Num label="L" value={selected.w} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { w: Math.max(0.1, v) })} />
          <Num label="B" value={selected.h} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { h: Math.max(0.1, v) })} />
          <Num label="H" value={selected.d} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { d: Math.max(0.1, v) })} />
        </Section>
      )}
      {selected?.kind === "circle" && (
        <Section title="Abmessungen">
          <Num label="cx" value={selected.cx} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { cx: v })} />
          <Num label="cy" value={selected.cy} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { cy: v })} />
          <Num label="r" value={selected.r} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { r: Math.max(0.1, v) })} />
        </Section>
      )}
      {selected?.kind === "cylinder" && (
        <Section title="Abmessungen">
          <Num label="cx" value={selected.cx} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { cx: v })} />
          <Num label="cy" value={selected.cy} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { cy: v })} />
          <Num label="cz" value={selected.cz} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { cz: v })} />
          <Num label="r" value={selected.r} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { r: Math.max(0.1, v) })} />
          <Num label="h" value={selected.height} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { height: Math.max(0.1, v) })} />
        </Section>
      )}
      {selected?.kind === "sphere" && (
        <Section title="Abmessungen">
          <Num label="cx" value={selected.cx} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { cx: v })} />
          <Num label="cy" value={selected.cy} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { cy: v })} />
          <Num label="cz" value={selected.cz} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { cz: v })} />
          <Num label="r" value={selected.r} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { r: Math.max(0.1, v) })} />
        </Section>
      )}
      <p className="text-xs text-subtle">
        Linke Maustaste dreht, mittleres schiebt, Rad zoomt. Zeichnen auf der XY-Ebene. F = einpassen.
      </p>
    </div>
  );
}

function MeshPanel() {
  const meshSize = usePre((s) => s.meshSize);
  const mesh = usePre((s) => s.mesh);
  const dim = usePre((s) => s.dim);
  const q = mesh ? meshQuality(mesh) : null;
  return (
    <div className="flex flex-col gap-3">
      <Section title="Vernetzung">
        <Num label="Kantenlänge" value={meshSize} unit="mm" step={0.5} onChange={(v) => usePre.getState().setMeshSize(Math.max(0.5, v))} />
        <p className="text-xs text-subtle">{dim === "3d" ? "C3D8 Hexeder, C3D6 Prismen." : "CPS4 Quadrate, CPS3 Delaunay."}</p>
        <Button className="mt-2 w-full" onClick={() => usePre.getState().generateMesh()}>
          <Layers />
          Netz erzeugen
        </Button>
      </Section>
      {q && (
        <Section title="Qualität">
          <Row k="Knoten" v={String(q.nnode)} />
          <Row k="Elemente" v={String(q.nelem)} />
          <Row k="Typen" v={q.types.join(", ")} />
          <Row k="min. Winkel" v={`${q.minAngle.toFixed(1)}°`} />
          <Row k="max. Aspekt" v={q.maxAspect.toFixed(2)} />
          <Row k="kritisch" v={String(q.nBad)} />
        </Section>
      )}
    </div>
  );
}

function MatPanel() {
  const materials = usePre((s) => s.materials);
  const selected = usePre((s) => s.shapes.find((x) => x.id === s.selectedShapeId));
  const dim = usePre((s) => s.dim);
  return (
    <div className="flex flex-col gap-3">
      <Section title="Zuweisung">
        {!selected && <p className="text-xs text-muted">Körper wählen, dann Material setzen.</p>}
        {selected && (
          <select
            value={selected.materialId}
            onChange={(e) => usePre.getState().setShapeMaterial(selected.id, e.target.value)}
            className="h-9 w-full rounded-md border border-border bg-surface-2 px-2 text-xs"
          >
            {materials.map((m) => (
              <option key={m.id} value={m.id}>
                {m.name}
              </option>
            ))}
          </select>
        )}
      </Section>
      {materials.map((m) => (
        <Section key={m.id} title={m.name}>
          <Num label="E" value={m.E} unit="MPa" onChange={(v) => usePre.getState().updateMaterial(m.id, { E: v })} />
          <Num label="ν" value={m.nu} unit="" step={0.01} onChange={(v) => usePre.getState().updateMaterial(m.id, { nu: v })} />
          <Num label="ρ" value={m.density} unit="t/mm³" step={1e-10} onChange={(v) => usePre.getState().updateMaterial(m.id, { density: v })} />
          {dim === "2d" && (
            <Num label="Dicke" value={m.thickness} unit="mm" step={0.1} onChange={(v) => usePre.getState().updateMaterial(m.id, { thickness: v })} />
          )}
        </Section>
      ))}
      <Button variant="secondary" size="sm" onClick={() => usePre.getState().addCustomMaterial()}>
        <Plus />
        Werkstoff
      </Button>
    </div>
  );
}

function BcPanel() {
  const dim = usePre((s) => s.dim);
  const restraints = usePre((s) => s.restraints);
  const loads = usePre((s) => s.loads);
  const selectedNodeIds = usePre((s) => s.selectedNodeIds);
  const selectedFace = usePre((s) => s.selectedFace);
  const gravity = loads.some((l) => l.kind === "gravity");
  const [fx, setFx] = useState(0);
  const [fy, setFy] = useState(0);
  const [fz, setFz] = useState(-200);
  const solid = dim === "3d";
  return (
    <div className="flex flex-col gap-3">
      <Section title={solid ? "Lager an Fläche" : "Lager an Kante"}>
        <div className="flex flex-wrap gap-1.5">
          {FACES.filter((f) => solid || f.id === "xmin" || f.id === "xmax" || f.id === "ymin" || f.id === "ymax").map((f) => (
            <Button
              key={f.id}
              variant="secondary"
              size="sm"
              onClick={() =>
                solid
                  ? usePre.getState().addFaceRestraint(f.id, true, true, true)
                  : usePre.getState().addEdgeRestraint(f.id === "xmin" ? "left" : f.id === "xmax" ? "right" : f.id === "ymin" ? "bottom" : "top", true, true)
              }
            >
              {f.label} fest
            </Button>
          ))}
        </div>
        {selectedFace && solid && (
          <Button className="mt-2 w-full" variant="secondary" size="sm" onClick={() => usePre.getState().addFaceRestraint(selectedFace.face, true, true, true)}>
            Fläche {faceLabel(selectedFace.face)} fest
          </Button>
        )}
        {selectedNodeIds.length > 0 && (
          <Button className="mt-2 w-full" variant="secondary" size="sm" onClick={() => usePre.getState().addNodeRestraint(true, true, solid)}>
            {selectedNodeIds.length} Knoten fest
          </Button>
        )}
      </Section>
      <Section title="Lager">
        {restraints.length === 0 && <p className="text-xs text-muted">Noch keine Lager.</p>}
        {restraints.map((r) => (
          <div key={r.id} className="flex items-center justify-between gap-2 text-xs">
            <span className="text-muted">
              {r.target.type === "face"
                ? faceLabel(r.target.face)
                : r.target.type === "edge"
                  ? faceLabel(r.target.edge)
                  : `${r.target.nodeIds.length} Knoten`}{" "}
              · {r.ux ? "Ux " : ""}
              {r.uy ? "Uy " : ""}
              {r.uz ? "Uz" : ""}
            </span>
            <button type="button" className="text-subtle hover:text-fg" onClick={() => usePre.getState().removeRestraint(r.id)}>
              ×
            </button>
          </div>
        ))}
      </Section>
      <Section title="Kraft (gesamt)">
        <Num label="Fx" value={fx} unit="N" onChange={setFx} />
        <Num label="Fy" value={fy} unit="N" onChange={setFy} />
        {solid && <Num label="Fz" value={fz} unit="N" onChange={setFz} />}
        <div className="mt-2 flex flex-wrap gap-1.5">
          {FACES.filter((f) => solid || f.id === "xmin" || f.id === "xmax" || f.id === "ymin" || f.id === "ymax").map((f) => (
            <Button
              key={f.id}
              variant="secondary"
              size="sm"
              onClick={() =>
                solid
                  ? usePre.getState().addFaceForce(f.id, fx, fy, fz)
                  : usePre.getState().addEdgeForce(
                      f.id === "xmin" ? "left" : f.id === "xmax" ? "right" : f.id === "ymin" ? "bottom" : "top",
                      fx,
                      fy,
                      0,
                    )
              }
            >
              {f.label}
            </Button>
          ))}
          {selectedNodeIds.length > 0 && (
            <Button variant="secondary" size="sm" onClick={() => usePre.getState().addNodeForce(fx, fy, solid ? fz : 0)}>
              Knoten
            </Button>
          )}
        </div>
      </Section>
      <Section title="Lasten">
        <label className="flex items-center gap-2 text-xs text-muted">
          <input type="checkbox" checked={gravity} onChange={(e) => usePre.getState().toggleGravity(e.target.checked)} className="size-4 accent-primary" />
          Eigengewicht {solid ? "(−Z)" : "(−Y)"}
        </label>
        {loads
          .filter((l) => l.kind === "force")
          .map((l) => (
            <div key={l.id} className="flex items-center justify-between gap-2 text-xs">
              <span className="text-muted">
                {l.kind === "force" && l.target.type === "face"
                  ? faceLabel(l.target.face)
                  : l.kind === "force" && l.target.type === "edge"
                    ? faceLabel(l.target.edge)
                    : "Knoten"}{" "}
                · {l.kind === "force" ? `${l.fx},${l.fy}${solid ? `,${l.fz}` : ""} N` : ""}
              </span>
              <button type="button" className="text-subtle hover:text-fg" onClick={() => usePre.getState().removeLoad(l.id)}>
                ×
              </button>
            </div>
          ))}
      </Section>
      <p className="text-xs text-subtle">Körper anklicken für Fläche. Werkzeug „Knoten“: Knoten wählen, dann Lager oder Kraft.</p>
    </div>
  );
}

function CheckPanel() {
  const issues = usePre((s) => s.issues);
  return (
    <div className="flex flex-col gap-3">
      <Section title="Modellcheck">
        <ul className="flex flex-col gap-1.5">
          {issues.map((i, idx) => (
            <li key={idx} className="flex gap-2 text-xs">
              <CheckCircle2
                className={cn("mt-0.5 size-3.5 shrink-0", i.level === "error" ? "text-danger" : i.level === "warn" ? "text-muted" : "text-ok")}
              />
              <span className={i.level === "error" ? "text-danger" : "text-muted"}>{i.message}</span>
            </li>
          ))}
        </ul>
      </Section>
    </div>
  );
}

function Section({ title, children }: { title: string; children: ReactNode }) {
  return (
    <section className="flex flex-col gap-1.5">
      <h2 className="text-xs font-medium tracking-wide text-muted uppercase">{title}</h2>
      {children}
    </section>
  );
}

function Num({
  label,
  value,
  onChange,
  unit,
  step = 1,
}: {
  label: string;
  value: number;
  onChange: (v: number) => void;
  unit: string;
  step?: number;
}) {
  return (
    <label className="flex items-center justify-between gap-2 py-0.5 text-xs">
      <span className="text-muted">{label}</span>
      <span className="flex items-center gap-1">
        <input
          type="number"
          step={step}
          value={Number.isFinite(value) ? value : 0}
          onChange={(e) => onChange(e.target.valueAsNumber)}
          className="h-8 w-24 rounded-md border border-border bg-surface-2 px-2 text-right font-mono text-xs"
        />
        <span className="w-10 text-subtle">{unit}</span>
      </span>
    </label>
  );
}

function Row({ k, v }: { k: string; v: string }) {
  return (
    <div className="flex justify-between gap-2 font-mono text-xs">
      <span className="text-muted">{k}</span>
      <span>{v}</span>
    </div>
  );
}
