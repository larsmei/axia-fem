/** Sequential cool-to-warm field map (steel → rust). */
const STOPS: [number, number, number][] = [
  [0.082, 0.141, 0.196],
  [0.165, 0.353, 0.447],
  [0.478, 0.62, 0.655],
  [0.82, 0.69, 0.545],
  [0.769, 0.361, 0.243],
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
