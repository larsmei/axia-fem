import { createFileRoute } from "@tanstack/react-router";
import { useCallback, useEffect, useMemo, useRef, useState, type ReactNode } from "react";
import {
  Download,
  FolderOpen,
  Play,
  RotateCcw,
} from "lucide-react";
import { Group, Panel, Separator } from "react-resizable-panels";
import { Button } from "@/components/ui/button";
import { InpEditor } from "@/components/inp-editor";
import { MeshViewer, fieldRange, type FieldId } from "@/components/mesh-viewer";
import { COLORBAR_CSS } from "@/lib/colormap";
import { DEFAULT_EXAMPLE, EXAMPLES } from "@/lib/examples";
import {
  initSolver,
  previewInp,
  solveInp,
  type FemResult,
} from "@/lib/solver";
import { cn, downloadText, formatNum } from "@/lib/utils";

export const Route = createFileRoute("/")({ component: Home });

const FIELDS: { id: FieldId; label: string }[] = [
  { id: "vm", label: "von Mises" },
  { id: "u", label: "|u|" },
  { id: "ux", label: "ux" },
  { id: "uy", label: "uy" },
  { id: "uz", label: "uz" },
  { id: "sxx", label: "Sxx" },
  { id: "syy", label: "Syy" },
  { id: "szz", label: "Szz" },
  { id: "sxy", label: "Sxy" },
];

function Home() {
  const [inp, setInp] = useState(DEFAULT_EXAMPLE.inp);
  const [exampleId, setExampleId] = useState(DEFAULT_EXAMPLE.id);
  const [ready, setReady] = useState(false);
  const [busy, setBusy] = useState(false);
  const [mesh, setMesh] = useState<FemResult | null>(null);
  const [result, setResult] = useState<FemResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [field, setField] = useState<FieldId>("vm");
  const [deformed, setDeformed] = useState(true);
  const [scaleMul, setScaleMul] = useState(1);
  const [tab, setTab] = useState<"inp" | "view">("view");
  const [log, setLog] = useState<string[]>(["Axia · FEM · C3D8 / T3D2 / S4R / B32 · NLGEOM / *PLASTIC"]);
  const fileRef = useRef<HTMLInputElement>(null);
  const autoRan = useRef(false);

  useEffect(() => {
    let live = true;
    initSolver()
      .then(() => {
        if (live) setReady(true);
      })
      .catch((e: unknown) => {
        if (live) setError(e instanceof Error ? e.message : "WASM konnte nicht geladen werden.");
      });
    return () => {
      live = false;
    };
  }, []);

  useEffect(() => {
    if (!ready || autoRan.current) return;
    autoRan.current = true;
    const r = solveInp(inp);
    if (!r.ok) {
      setError(r.error ?? "Rechnung fehlgeschlagen.");
      return;
    }
    setResult(r);
    setMesh(r);
    const s = r.stats;
    if (s) {
      setLog((l) => [
        `${s.procedure ? s.procedure + " · " : ""}${s.solver} · ${s.iterations} it · ${s.nnode} Knoten · ${s.nelem} Elemente · ${s.nfree}/${s.ndof} DOF · |u|max ${formatNum(s.uMax, 5)} · σvm ${formatNum(s.vmMax, 3)}${s.lambda != null && s.procedure?.includes("RIKS") ? ` · λ ${formatNum(s.lambda, 4)}` : ""} · ${s.timeMs.toFixed(0)} ms`,
        ...l,
      ].slice(0, 12));
    }
  }, [ready, inp]);

  useEffect(() => {
    if (!ready) return;
    const t = window.setTimeout(() => {
      const p = previewInp(inp);
      if (p.ok) {
        setMesh(p);
        setError(null);
      } else {
        setMesh(p);
        setError(p.error ?? "INP nicht lesbar.");
      }
    }, 280);
    return () => window.clearTimeout(t);
  }, [inp, ready]);

  const run = useCallback(() => {
    if (!ready || busy) return;
    setBusy(true);
    setError(null);
    window.setTimeout(() => {
      const r = solveInp(inp);
      setBusy(false);
      if (!r.ok) {
        setError(r.error ?? "Rechnung fehlgeschlagen.");
        setLog((l) => [r.error ?? "Fehler", ...l].slice(0, 12));
        return;
      }
      setResult(r);
      setMesh(r);
      const s = r.stats;
      const line = s
        ? `${s.procedure ? s.procedure + " · " : ""}${s.solver} · ${s.iterations} it · ${s.nnode} Knoten · ${s.nelem} Elemente · ${s.nfree}/${s.ndof} DOF · |u|max ${formatNum(s.uMax, 5)} · σvm ${formatNum(s.vmMax, 3)}${s.lambda != null && s.procedure?.includes("RIKS") ? ` · λ ${formatNum(s.lambda, 4)}` : ""} · ${s.timeMs.toFixed(0)} ms`
        : "Fertig.";
      setLog((l) => [line, ...(r.warnings ?? []), ...l].slice(0, 12));
      setTab("view");
    }, 30);
  }, [busy, inp, ready]);

  const loadExample = (id: string) => {
    const ex = EXAMPLES.find((e) => e.id === id);
    if (!ex) return;
    setExampleId(id);
    setInp(ex.inp);
    setError(null);
    if (!ready) {
      setResult(null);
      return;
    }
    const r = solveInp(ex.inp);
    if (!r.ok) {
      setResult(null);
      setError(r.error ?? "Rechnung fehlgeschlagen.");
      return;
    }
    setResult(r);
    setMesh(r);
    const s = r.stats;
    if (s) {
      setLog((l) =>
        [
          `${s.procedure ? s.procedure + " · " : ""}${s.solver} · ${s.iterations} it · ${s.nnode} Knoten · ${s.nelem} Elemente · ${s.nfree}/${s.ndof} DOF · |u|max ${formatNum(s.uMax, 5)} · σvm ${formatNum(s.vmMax, 3)}${s.lambda != null && s.procedure?.includes("RIKS") ? ` · λ ${formatNum(s.lambda, 4)}` : ""} · ${s.timeMs.toFixed(0)} ms`,
          ...l,
        ].slice(0, 12),
      );
    }
  };

  const onFile = (file: File) => {
    file.text().then((t) => {
      setExampleId("file");
      setInp(t);
      setResult(null);
    });
  };

  const autoScale = useMemo(() => {
    const um = result?.stats?.uMax ?? 0;
    if (!um || !mesh?.coords?.length) return 1;
    let span = 0;
    const c = mesh.coords;
    let minx = Infinity,
      maxx = -Infinity,
      miny = Infinity,
      maxy = -Infinity,
      minz = Infinity,
      maxz = -Infinity;
    for (let i = 0; i < c.length; i += 3) {
      minx = Math.min(minx, c[i]);
      maxx = Math.max(maxx, c[i]);
      miny = Math.min(miny, c[i + 1]);
      maxy = Math.max(maxy, c[i + 1]);
      minz = Math.min(minz, c[i + 2]);
      maxz = Math.max(maxz, c[i + 2]);
    }
    span = Math.max(maxx - minx, maxy - miny, maxz - minz, 1);
    return (0.14 * span) / um;
  }, [mesh, result]);

  const scale = autoScale * scaleMul;
  const range = fieldRange(mesh, result, field);
  const solved = Boolean(result?.ok && result.kind === "solve");

  return (
    <main className="flex h-dvh flex-col overflow-hidden bg-bg text-fg">
      <header className="flex shrink-0 items-center gap-2 border-b border-border bg-surface px-3 py-2 lg:px-4">
        <Mark />
        <div className="min-w-0 flex-1">
          <div className="flex items-baseline gap-2">
            <h1 className="text-[15px] font-medium tracking-tight">Axia</h1>
            <p className="hidden text-xs text-muted sm:block">FEM-Solver · INP / FRD</p>
          </div>
        </div>
        <label className="hidden sm:block">
          <span className="sr-only">Beispiel</span>
          <select
            value={exampleId}
            onChange={(e) => loadExample(e.target.value)}
            className="h-10 max-w-[11rem] rounded-md border border-border bg-surface-2 px-2 text-sm text-fg"
          >
            {EXAMPLES.map((e) => (
              <option key={e.id} value={e.id}>
                {e.name}
              </option>
            ))}
            {exampleId === "file" && <option value="file">Datei</option>}
          </select>
        </label>
        <input
          ref={fileRef}
          type="file"
          accept=".inp,.txt"
          className="hidden"
          onChange={(e) => {
            const f = e.target.files?.[0];
            if (f) onFile(f);
            e.target.value = "";
          }}
        />
        <Button variant="secondary" size="icon" className="lg:hidden" onClick={() => fileRef.current?.click()} aria-label="INP öffnen">
          <FolderOpen />
        </Button>
        <Button variant="secondary" size="sm" className="hidden lg:inline-flex" onClick={() => fileRef.current?.click()}>
          <FolderOpen />
          Öffnen
        </Button>
        <Button onClick={run} disabled={!ready || busy} className="min-w-24">
          <Play />
          {busy ? "Rechnet…" : "Lösen"}
        </Button>
      </header>

      <div className="flex gap-1 border-b border-border px-2 py-1 lg:hidden">
        <TabBtn active={tab === "inp"} onClick={() => setTab("inp")}>
          INP
        </TabBtn>
        <TabBtn active={tab === "view"} onClick={() => setTab("view")}>
          Netz
        </TabBtn>
        <select
          value={exampleId}
          onChange={(e) => loadExample(e.target.value)}
          className="ml-auto h-10 max-w-[10rem] rounded-md border border-border bg-surface-2 px-2 text-sm"
        >
          {EXAMPLES.map((e) => (
            <option key={e.id} value={e.id}>
              {e.name}
            </option>
          ))}
        </select>
      </div>

      <div className="flex min-h-0 flex-1">
        <div className="hidden min-h-0 flex-1 lg:block">
          <Group orientation="horizontal" className="h-full">
            <Panel defaultSize="38" minSize="22" className="min-h-0">
              <EditorColumn
                inp={inp}
                setInp={(v) => {
                  setInp(v);
                  setResult(null);
                }}
                log={log}
                error={error}
                mesh={mesh}
                result={result}
                ready={ready}
              />
            </Panel>
            <Separator className="w-px bg-border hover:bg-fg/30" />
            <Panel defaultSize="62" minSize="30" className="min-h-0">
              <ViewerColumn
                mesh={mesh}
                result={result}
                field={field}
                setField={setField}
                deformed={deformed}
                setDeformed={setDeformed}
                scaleMul={scaleMul}
                setScaleMul={setScaleMul}
                scale={scale}
                range={range}
                solved={solved}
                error={error}
              />
            </Panel>
          </Group>
        </div>

        <div className="flex min-h-0 min-w-0 flex-1 flex-col lg:hidden">
          {tab === "inp" ? (
            <EditorColumn
              inp={inp}
              setInp={(v) => {
                setInp(v);
                setResult(null);
              }}
              log={log}
              error={error}
              mesh={mesh}
              result={result}
              ready={ready}
            />
          ) : (
            <ViewerColumn
              mesh={mesh}
              result={result}
              field={field}
              setField={setField}
              deformed={deformed}
              setDeformed={setDeformed}
              scaleMul={scaleMul}
              setScaleMul={setScaleMul}
              scale={scale}
              range={range}
              solved={solved}
              error={error}
            />
          )}
        </div>
      </div>
    </main>
  );
}

function Mark() {
  return (
    <svg width="28" height="28" viewBox="0 0 28 28" aria-hidden className="shrink-0">
      <rect width="28" height="28" rx="7" className="fill-surface-2 stroke-border" strokeWidth="1" />
      <path
        d="M8 20 L14 7 L20 20 Z"
        fill="none"
        className="stroke-fg"
        strokeWidth="1.4"
        strokeLinejoin="round"
      />
      <circle cx="14" cy="7" r="1.4" className="fill-fg" />
      <circle cx="8" cy="20" r="1.4" className="fill-fg" />
      <circle cx="20" cy="20" r="1.4" className="fill-fg" />
    </svg>
  );
}

function TabBtn({
  active,
  onClick,
  children,
}: {
  active: boolean;
  onClick: () => void;
  children: ReactNode;
}) {
  return (
    <button
      type="button"
      onClick={onClick}
      className={cn(
        "h-10 min-w-16 rounded-md px-3 text-sm font-medium",
        active ? "bg-surface-2 text-fg" : "text-muted",
      )}
    >
      {children}
    </button>
  );
}

function EditorColumn({
  inp,
  setInp,
  log,
  error,
  mesh,
  result,
  ready,
}: {
  inp: string;
  setInp: (v: string) => void;
  log: string[];
  error: string | null;
  mesh: FemResult | null;
  result: FemResult | null;
  ready: boolean;
}) {
  return (
    <div className="flex h-full min-h-0 flex-col bg-surface">
      <div className="flex items-center justify-between border-b border-border px-3 py-2">
        <p className="text-xs font-medium tracking-wide text-muted uppercase">Eingabe · job.inp</p>
        <p className="font-mono text-[11px] text-subtle tabular-nums">
          {ready ? `${mesh?.nnode ?? "—"} nd · ${mesh?.nelem ?? "—"} el` : "lädt Solver…"}
        </p>
      </div>
      <InpEditor value={inp} onChange={setInp} />
      <div className="flex shrink-0 gap-2 border-t border-border p-2">
        <Button
          variant="secondary"
          size="sm"
          disabled={!result?.frd}
          onClick={() => result?.frd && downloadText("job.frd", result.frd, "text/plain")}
        >
          <Download />
          FRD
        </Button>
        <Button
          variant="secondary"
          size="sm"
          disabled={!result?.dat}
          onClick={() => result?.dat && downloadText("job.dat", result.dat, "text/plain")}
        >
          <Download />
          DAT
        </Button>
        <Button variant="ghost" size="sm" onClick={() => downloadText("job.inp", inp, "text/plain")}>
          INP
        </Button>
      </div>
      <div className="h-[7.5rem] shrink-0 overflow-auto border-t border-border bg-code px-3 py-2 font-mono text-[11px] leading-5 text-muted">
        {error && <p className="text-danger">{error}</p>}
        {log.map((line, i) => (
          <p key={i} className={i === 0 ? "text-fg" : ""}>
            {line}
          </p>
        ))}
      </div>
    </div>
  );
}

function ViewerColumn({
  mesh,
  result,
  field,
  setField,
  deformed,
  setDeformed,
  scaleMul,
  setScaleMul,
  scale,
  range,
  solved,
  error,
}: {
  mesh: FemResult | null;
  result: FemResult | null;
  field: FieldId;
  setField: (f: FieldId) => void;
  deformed: boolean;
  setDeformed: (v: boolean) => void;
  scaleMul: number;
  setScaleMul: (v: number) => void;
  scale: number;
  range: { min: number; max: number } | null;
  solved: boolean;
  error: string | null;
}) {
  return (
    <div className="relative flex h-full min-h-0 flex-col bg-viewport">
      <div className="pointer-events-none absolute inset-x-0 top-0 z-10 flex flex-wrap items-start justify-between gap-2 p-3">
        <div className="pointer-events-auto flex flex-wrap items-center gap-1.5 rounded-xl bg-bg/80 p-1.5 shadow-[0_0_0_1px_rgba(255,255,255,0.08)]">
          {FIELDS.map((f) => (
            <button
              key={f.id}
              type="button"
              disabled={!solved}
              onClick={() => setField(f.id)}
              className={cn(
                "h-8 rounded-md px-2.5 text-xs font-medium",
                field === f.id && solved ? "bg-primary text-primary-foreground" : "text-muted hover:text-fg",
              )}
            >
              {f.label}
            </button>
          ))}
        </div>
        <div className="pointer-events-auto flex items-center gap-2 rounded-xl bg-bg/80 px-2.5 py-1.5 shadow-[0_0_0_1px_rgba(255,255,255,0.08)]">
          <label className="flex items-center gap-2 text-xs text-muted">
            <input
              type="checkbox"
              checked={deformed}
              onChange={(e) => setDeformed(e.target.checked)}
              className="size-4 accent-primary"
            />
            Verformt
          </label>
          <input
            type="range"
            min={0}
            max={30}
            step={0.1}
            value={scaleMul}
            onChange={(e) => setScaleMul(Number(e.target.value))}
            className="w-24"
            aria-label="Verformungsskalierung"
          />
          <span className="w-10 font-mono text-[11px] text-subtle tabular-nums">{formatNum(scale, 2)}×</span>
          <button
            type="button"
            className="text-muted hover:text-fg"
            onClick={() => setScaleMul(1)}
            aria-label="Skalierung zurücksetzen"
          >
            <RotateCcw className="size-3.5" />
          </button>
        </div>
      </div>

      <MeshViewer mesh={mesh} result={solved ? result : mesh} field={field} deformed={deformed && solved} scale={scale} />

      <div className="absolute right-3 bottom-16 z-10 flex items-stretch gap-2">
        <div
          className="h-36 w-2.5 rounded-full"
          style={{ background: `linear-gradient(to top, ${COLORBAR_CSS})` }}
        />
        <div className="flex flex-col justify-between py-0.5 font-mono text-[11px] text-muted tabular-nums">
          <span>{solved && range ? formatNum(range.max, 3) : "—"}</span>
          <span>{solved && range ? formatNum((range.min + range.max) / 2, 3) : ""}</span>
          <span>{solved && range ? formatNum(range.min, 3) : "—"}</span>
        </div>
      </div>

      <footer className="flex shrink-0 items-center justify-between gap-3 border-t border-border bg-surface/90 px-3 py-2 text-xs text-muted">
        <p className="min-w-0 truncate">
          {error
            ? error
            : solved && result?.stats
              ? `${result.stats.solver} · u_max ${formatNum(result.stats.uMax, 5)} · σ_vm ${formatNum(result.stats.vmMin, 2)} … ${formatNum(result.stats.vmMax, 2)}`
              : mesh?.ok
                ? `${mesh.nnode} Knoten, ${mesh.nelem} Elemente — Lösen startet die Rechnung.`
                : "INP wird gelesen."}
        </p>
        <p className="hidden shrink-0 sm:block">{mesh?.heading}</p>
      </footer>
    </div>
  );
}
