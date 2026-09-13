/** Classic jet / rainbow (Mecway, MATLAB, cgx). Low = navy, high = dark red. */
const STOPS: [number, number, number][] = [
  [0.0, 0.0, 0.5],
  [0.0, 0.0, 1.0],
  [0.0, 0.5, 1.0],
  [0.0, 1.0, 1.0],
  [0.5, 1.0, 0.5],
  [1.0, 1.0, 0.0],
  [1.0, 0.5, 0.0],
  [1.0, 0.0, 0.0],
  [0.5, 0.0, 0.0],
];

function lerp(a: number, b: number, t: number) {
  return a + (b - a) * t;
}

export function sampleColor(t: number): [number, number, number] {
  const x = Math.min(1, Math.max(0, t));
  const n = STOPS.length - 1;
  const p = x * n;
  const i = Math.min(n - 1, Math.floor(p));
  const f = p - i;
  const a = STOPS[i];
  const b = STOPS[i + 1];
  return [lerp(a[0], b[0], f), lerp(a[1], b[1], f), lerp(a[2], b[2], f)];
}

export function sampleCss(t: number): string {
  const [r, g, b] = sampleColor(t);
  return `rgb(${Math.round(r * 255)}, ${Math.round(g * 255)}, ${Math.round(b * 255)})`;
}

export const COLORBAR_CSS = STOPS.map(
  ([r, g, b]) => `rgb(${Math.round(r * 255)} ${Math.round(g * 255)} ${Math.round(b * 255)})`,
).join(", ");

export const COLORBAR_TICKS = 10;

/** Mecway-style scientific ticks: `1.122E+07`, small integers plain. */
export function formatLegend(v: number): string {
  if (!Number.isFinite(v)) return "—";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 100 && a < 10000) {
    const r = Math.round(v);
    if (Math.abs(v - r) <= Math.max(a * 1e-4, 0.51)) return String(r);
  }
  const [mant, expRaw] = v.toExponential(3).split("e");
  const e = Number(expRaw);
  const sign = e >= 0 ? "+" : "-";
  return `${mant}E${sign}${String(Math.abs(e)).padStart(2, "0")}`;
}
