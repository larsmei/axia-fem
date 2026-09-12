export type FemElement = {
  id: number;
  type: string;
  nodes: number[];
  secA?: number;
  secB?: number;
  n1?: number[];
  area?: number;
  th?: number;
};

export type FemStats = {
  nnode: number;
  nelem: number;
  ndof: number;
  nfree: number;
  solver: string;
  iterations: number;
  residual: number;
  timeMs: number;
  uMax: number;
  vmMin: number;
  vmMax: number;
  nbc: number;
  ncload: number;
};

export type FemResult = {
  ok: boolean;
  kind: "preview" | "solve" | "error";
  error?: string;
  heading?: string;
  dim?: number;
  nodeIds?: number[];
  coords?: number[];
  elements?: FemElement[];
  materials?: { name: string; E: number; nu: number; density: number }[];
  nnode?: number;
  nelem?: number;
  warnings?: string[];
  u?: number[];
  ur?: number[];
  stress?: number[];
  rf?: number[];
  vonMises?: number[];
  frd?: string;
  dat?: string;
  stats?: FemStats;
};

type WasmMod = {
  default: (opts?: { module_or_path?: string }) => Promise<unknown>;
  preview_inp: (inp: string) => string;
  solve_inp: (inp: string) => string;
};

let api: WasmMod | null = null;
let loading: Promise<void> | null = null;

export function isSolverReady() {
  return api !== null;
}

export async function initSolver() {
  if (api) return;
  if (!loading) {
    loading = (async () => {
      const [mod] = await Promise.all([
        import("@/wasm/axia_fem.js") as Promise<WasmMod>,
      ]);
      await mod.default({ module_or_path: "/axia_fem_bg.wasm" });
      api = mod;
    })();
  }
  await loading;
}

function parse(json: string): FemResult {
  try {
    return JSON.parse(json) as FemResult;
  } catch {
    return { ok: false, kind: "error", error: "Ungültige Solver-Antwort." };
  }
}

export function previewInp(inp: string): FemResult {
  if (!api) return { ok: false, kind: "error", error: "Solver nicht geladen." };
  return parse(api.preview_inp(inp));
}

export function solveInp(inp: string): FemResult {
  if (!api) return { ok: false, kind: "error", error: "Solver nicht geladen." };
  const t0 = performance.now();
  const r = parse(api.solve_inp(inp));
  if (r.ok && r.stats) {
    r.stats = { ...r.stats, timeMs: performance.now() - t0 };
  }
  return r;
}
