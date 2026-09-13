import { useEffect, useMemo, useState, type ReactNode } from "react";
import { useNavigate } from "@tanstack/react-router";
import {
  BoxSelect,
  CheckCircle2,
  Circle,
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
import { meshQuality } from "@/lib/pre/mesh";
import { usePre } from "@/lib/pre/store";
import { TEMPLATE_META } from "@/lib/pre/templates";
import type { EdgeName, Tool } from "@/lib/pre/types";
import { cn, downloadText } from "@/lib/utils";

type Step = "geo" | "mesh" | "mat" | "bc" | "check";

const STEPS: { id: Step; n: string; label: string }[] = [
  { id: "geo", n: "1", label: "Geometrie" },
  { id: "mesh", n: "2", label: "Netz" },
  { id: "mat", n: "3", label: "Material" },
  { id: "bc", n: "4", label: "Lager / Last" },
  { id: "check", n: "5", label: "Prüfen" },
];

export function Preprocessor() {
  const navigate = useNavigate();
  const [step, setStep] = useState<Step>("geo");
  const name = usePre((s) => s.name);
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
            <p className="hidden text-xs text-muted sm:block">Präprozessor · 2D CPS</p>
          </div>
        </div>
        <AppNav active="pre" />
        <input
          value={name}
          onChange={(e) => usePre.getState().setName(e.target.value)}
          className="hidden h-8 w-36 rounded-md border border-border bg-surface-2 px-2 text-xs text-fg sm:block"
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
          <ToolBtn id="rect" icon={Square} label="Rechteck" />
          <ToolBtn id="circle" icon={Circle} label="Kreis" />
          <ToolBtn id="polygon" icon={Hexagon} label="Polygon" />
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
            <ToolBtn id="rect" icon={Square} label="Rechteck" wide />
            <ToolBtn id="circle" icon={Circle} label="Kreis" wide />
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
          {shapes.length} Körper · {q ? `${q.nnode} Knoten · ${q.nelem} Elemente` : "kein Netz"} · mm / N / MPa
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
  const shapes = usePre((s) => s.shapes);
  const selected = usePre((s) => s.shapes.find((x) => x.id === s.selectedShapeId));
  return (
    <div className="flex flex-col gap-3">
      <Section title="Vorlagen">
        <div className="flex flex-wrap gap-1.5">
          {TEMPLATE_META.map((t) => (
            <Button key={t.id} variant="secondary" size="sm" onClick={() => usePre.getState().loadTemplate(t.id)}>
              {t.name}
            </Button>
          ))}
          <Button variant="ghost" size="sm" onClick={() => usePre.getState().reset()}>
            Neu
          </Button>
        </div>
      </Section>
      <Section title="Körper">
        {shapes.length === 0 && <p className="text-xs text-muted">Rechteck, Kreis oder Polygon zeichnen.</p>}
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
                {s.kind === "rect" ? "Rechteck" : s.kind === "circle" ? "Kreis" : "Polygon"} {i + 1}
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
        </Section>
      )}
      {selected?.kind === "circle" && (
        <Section title="Abmessungen">
          <Num label="cx" value={selected.cx} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { cx: v })} />
          <Num label="cy" value={selected.cy} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { cy: v })} />
          <Num label="r" value={selected.r} unit="mm" onChange={(v) => usePre.getState().updateShape(selected.id, { r: Math.max(0.1, v) })} />
        </Section>
      )}
      <p className="text-xs text-subtle">Mausrad zoomt, mittleres Klicken schiebt. F = einpassen. Esc bricht ab.</p>
    </div>
  );
}

function MeshPanel() {
  const meshSize = usePre((s) => s.meshSize);
  const mesh = usePre((s) => s.mesh);
  const q = mesh ? meshQuality(mesh) : null;
  return (
    <div className="flex flex-col gap-3">
      <Section title="Vernetzung">
        <Num label="Kantenlänge" value={meshSize} unit="mm" step={0.5} onChange={(v) => usePre.getState().setMeshSize(Math.max(0.5, v))} />
        <Button className="mt-2 w-full" onClick={() => usePre.getState().generateMesh()}>
          <Layers />
          Netz erzeugen
        </Button>
      </Section>
      {q && (
        <Section title="Qualität">
          <Row k="Knoten" v={String(q.nnode)} />
          <Row k="Elemente" v={String(q.nelem)} />
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
          <Num label="Dicke" value={m.thickness} unit="mm" step={0.1} onChange={(v) => usePre.getState().updateMaterial(m.id, { thickness: v })} />
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
  const restraints = usePre((s) => s.restraints);
  const loads = usePre((s) => s.loads);
  const selectedNodeIds = usePre((s) => s.selectedNodeIds);
  const gravity = loads.some((l) => l.kind === "gravity");
  const [fx, setFx] = useState(0);
  const [fy, setFy] = useState(-100);
  const edge = (e: EdgeName, label: string) => (
    <Button key={e} variant="secondary" size="sm" onClick={() => usePre.getState().addEdgeRestraint(e, true, true)}>
      {label}
    </Button>
  );
  return (
    <div className="flex flex-col gap-3">
      <Section title="Lager an Kante">
        <div className="flex flex-wrap gap-1.5">
          {edge("left", "Links fest")}
          {edge("right", "Rechts fest")}
          {edge("bottom", "Unten fest")}
          {edge("top", "Oben fest")}
        </div>
        {selectedNodeIds.length > 0 && (
          <Button className="mt-2 w-full" variant="secondary" size="sm" onClick={() => usePre.getState().addNodeRestraint(true, true)}>
            {selectedNodeIds.length} Knoten fest
          </Button>
        )}
      </Section>
      <Section title="Lager">
        {restraints.length === 0 && <p className="text-xs text-muted">Noch keine Lager.</p>}
        {restraints.map((r) => (
          <div key={r.id} className="flex items-center justify-between gap-2 text-xs">
            <span className="text-muted">
              {r.target.type === "edge" ? r.target.edge : `${r.target.nodeIds.length} Knoten`} · {r.ux ? "Ux " : ""}
              {r.uy ? "Uy" : ""}
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
        <div className="mt-2 flex flex-wrap gap-1.5">
          <Button variant="secondary" size="sm" onClick={() => usePre.getState().addEdgeForce("right", fx, fy)}>
            Rechts
          </Button>
          <Button variant="secondary" size="sm" onClick={() => usePre.getState().addEdgeForce("top", fx, fy)}>
            Oben
          </Button>
          <Button variant="secondary" size="sm" onClick={() => usePre.getState().addEdgeForce("left", fx, fy)}>
            Links
          </Button>
          {selectedNodeIds.length > 0 && (
            <Button variant="secondary" size="sm" onClick={() => usePre.getState().addNodeForce(fx, fy)}>
              Knoten
            </Button>
          )}
        </div>
      </Section>
      <Section title="Lasten">
        <label className="flex items-center gap-2 text-xs text-muted">
          <input type="checkbox" checked={gravity} onChange={(e) => usePre.getState().toggleGravity(e.target.checked)} className="size-4 accent-primary" />
          Eigengewicht (g = 9810 mm/s²)
        </label>
        {loads
          .filter((l) => l.kind === "force")
          .map((l) => (
            <div key={l.id} className="flex items-center justify-between gap-2 text-xs">
              <span className="text-muted">
                {l.kind === "force" && l.target.type === "edge" ? l.target.edge : "Knoten"} · {l.kind === "force" ? `${l.fx},${l.fy} N` : ""}
              </span>
              <button type="button" className="text-subtle hover:text-fg" onClick={() => usePre.getState().removeLoad(l.id)}>
                ×
              </button>
            </div>
          ))}
      </Section>
      <p className="text-xs text-subtle">Werkzeug „Knoten“: Knoten klicken, dann Lager oder Kraft zuweisen.</p>
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
