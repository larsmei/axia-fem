import { useCallback, useEffect, useRef } from "react";
import { resolveNodes } from "@/lib/pre/export-inp";
import { dist, hitShape, modelBounds, nearestNode, shapeHoles, shapeOutline, snapPt } from "@/lib/pre/geometry";
import { usePre } from "@/lib/pre/store";
import type { Vec2 } from "@/lib/pre/types";

const BG = "#f4f5f7";
const GRID = "#e4e4e7";
const MAJOR = "#d4d4d8";
const INK = "#18181c";
const MUTED = "#71717a";
const FILL = "rgba(24,24,28,0.07)";
const SEL = "#09090b";
const MESH = "#8b8b94";
const DANGER = "#d4786a";

type Cam = { x: number; y: number; k: number };

export function PreViewport() {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const wrapRef = useRef<HTMLDivElement>(null);
  const cam = useRef<Cam>({ x: 50, y: 20, k: 6 });
  const drag = useRef<
    | { mode: "pan"; lx: number; ly: number }
    | { mode: "draw"; }
    | { mode: "box"; ax: number; ay: number; bx: number; by: number }
    | null
  >(null);
  const shapes = usePre((s) => s.shapes);
  const mesh = usePre((s) => s.mesh);
  const tool = usePre((s) => s.tool);
  const draft = usePre((s) => s.draft);
  const selectedShapeId = usePre((s) => s.selectedShapeId);
  const selectedNodeIds = usePre((s) => s.selectedNodeIds);
  const restraints = usePre((s) => s.restraints);
  const loads = usePre((s) => s.loads);
  const issues = usePre((s) => s.issues);

  const world = useCallback((cx: number, cy: number, w: number, h: number): Vec2 => {
    const c = cam.current;
    return { x: (cx - w / 2) / c.k + c.x, y: (h / 2 - cy) / c.k + c.y };
  }, []);

  const screen = (p: Vec2, w: number, h: number) => {
    const c = cam.current;
    return { x: (p.x - c.x) * c.k + w / 2, y: h / 2 - (p.y - c.y) * c.k };
  };

  const fit = useCallback(() => {
    const el = wrapRef.current;
    if (!el) return;
    const w = el.clientWidth;
    const h = el.clientHeight;
    const b = modelBounds(usePre.getState().shapes);
    const k = Math.min(w / Math.max(b.w, 1), h / Math.max(b.h, 1)) * 0.86;
    cam.current = { x: b.x + b.w / 2, y: b.y + b.h / 2, k: Math.max(2, Math.min(48, k)) };
  }, []);

  useEffect(() => {
    fit();
  }, [fit, shapes.length]);

  const paint = useCallback(() => {
    const canvas = canvasRef.current;
    const wrap = wrapRef.current;
    if (!canvas || !wrap) return;
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = wrap.clientWidth;
    const h = wrap.clientHeight;
    if (canvas.width !== Math.floor(w * dpr) || canvas.height !== Math.floor(h * dpr)) {
      canvas.width = Math.floor(w * dpr);
      canvas.height = Math.floor(h * dpr);
      canvas.style.width = `${w}px`;
      canvas.style.height = `${h}px`;
    }
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.fillStyle = BG;
    ctx.fillRect(0, 0, w, h);

    const c = cam.current;
    const step = c.k > 14 ? 1 : c.k > 6 ? 5 : 10;
    const major = step * 10;
    const w0 = world(0, h, w, h);
    const w1 = world(w, 0, w, h);
    ctx.lineWidth = 1;
    const x0 = Math.floor(w0.x / step) * step;
    const y0 = Math.floor(w0.y / step) * step;
    for (let x = x0; x <= w1.x + step; x += step) {
      const s = screen({ x, y: 0 }, w, h);
      ctx.strokeStyle = Math.abs(x % major) < 1e-6 ? MAJOR : GRID;
      ctx.beginPath();
      ctx.moveTo(Math.round(s.x) + 0.5, 0);
      ctx.lineTo(Math.round(s.x) + 0.5, h);
      ctx.stroke();
    }
    for (let y = y0; y <= w1.y + step; y += step) {
      const s = screen({ x: 0, y }, w, h);
      ctx.strokeStyle = Math.abs(y % major) < 1e-6 ? MAJOR : GRID;
      ctx.beginPath();
      ctx.moveTo(0, Math.round(s.y) + 0.5);
      ctx.lineTo(w, Math.round(s.y) + 0.5);
      ctx.stroke();
    }

    const ox = screen({ x: 0, y: 0 }, w, h);
    ctx.strokeStyle = INK;
    ctx.lineWidth = 1.25;
    ctx.beginPath();
    ctx.moveTo(ox.x, ox.y);
    ctx.lineTo(ox.x + 28, ox.y);
    ctx.moveTo(ox.x, ox.y);
    ctx.lineTo(ox.x, ox.y - 28);
    ctx.stroke();
    ctx.fillStyle = MUTED;
    ctx.font = "11px IBM Plex Sans, sans-serif";
    ctx.fillText("x", ox.x + 32, ox.y + 4);
    ctx.fillText("y", ox.x - 4, ox.y - 32);

    const st = usePre.getState();
    for (const sh of st.shapes) {
      const ring = shapeOutline(sh);
      pathRing(ctx, ring, w, h, screen);
      ctx.fillStyle = sh.id === st.selectedShapeId ? "rgba(24,24,28,0.12)" : FILL;
      ctx.fill();
      for (const hole of shapeHoles(sh)) {
        ctx.beginPath();
        const n = 32;
        for (let i = 0; i <= n; i++) {
          const a = (i / n) * Math.PI * 2;
          const p = screen({ x: hole.cx + hole.r * Math.cos(a), y: hole.cy + hole.r * Math.sin(a) }, w, h);
          if (i === 0) ctx.moveTo(p.x, p.y);
          else ctx.lineTo(p.x, p.y);
        }
        ctx.fillStyle = BG;
        ctx.fill();
      }
    }

    if (st.mesh) {
      const byId = new Map(st.mesh.nodes.map((n) => [n.id, n]));
      ctx.strokeStyle = MESH;
      ctx.lineWidth = 0.8;
      ctx.beginPath();
      for (const el of st.mesh.elements) {
        const pts = el.nodes.map((id) => byId.get(id)).filter(Boolean);
        if (pts.length < 3) continue;
        const p0 = screen(pts[0]!, w, h);
        ctx.moveTo(p0.x, p0.y);
        for (let i = 1; i < pts.length; i++) {
          const p = screen(pts[i]!, w, h);
          ctx.lineTo(p.x, p.y);
        }
        ctx.closePath();
      }
      ctx.stroke();
      if (c.k > 10) {
        ctx.fillStyle = INK;
        for (const n of st.mesh.nodes) {
          const p = screen(n, w, h);
          ctx.beginPath();
          ctx.arc(p.x, p.y, st.selectedNodeIds.includes(n.id) ? 3.5 : 1.6, 0, Math.PI * 2);
          ctx.fill();
        }
      }
    }

    for (const sh of st.shapes) {
      const ring = shapeOutline(sh);
      pathRing(ctx, ring, w, h, screen);
      ctx.strokeStyle = sh.id === st.selectedShapeId ? SEL : INK;
      ctx.lineWidth = sh.id === st.selectedShapeId ? 2.2 : 1.2;
      ctx.stroke();
    }

    if (st.mesh) {
      ctx.fillStyle = INK;
      for (const r of st.restraints) {
        const ids = resolveNodes(st.mesh, st.shapes, r.target);
        for (const id of ids) {
          const n = st.mesh.nodes.find((nd) => nd.id === id);
          if (!n) continue;
          const p = screen(n, w, h);
          drawSupport(ctx, p.x, p.y, r.ux, r.uy);
        }
      }
      ctx.strokeStyle = INK;
      ctx.fillStyle = INK;
      for (const ld of st.loads) {
        if (ld.kind !== "force") continue;
        const ids = resolveNodes(st.mesh, st.shapes, ld.target);
        for (const id of ids) {
          const n = st.mesh.nodes.find((nd) => nd.id === id);
          if (!n) continue;
          const p = screen(n, w, h);
          const mag = Math.hypot(ld.fx, ld.fy) || 1;
          const len = 18;
          drawArrow(ctx, p.x, p.y, (ld.fx / mag) * len, -(ld.fy / mag) * len);
        }
      }
    }

    if (st.draft?.tool === "rect") {
      const a = screen(st.draft.a, w, h);
      const b = screen(st.draft.b, w, h);
      ctx.strokeStyle = SEL;
      ctx.setLineDash([4, 3]);
      ctx.strokeRect(Math.min(a.x, b.x), Math.min(a.y, b.y), Math.abs(b.x - a.x), Math.abs(b.y - a.y));
      ctx.setLineDash([]);
      dimLabel(ctx, st.draft.a, st.draft.b, w, h, screen);
    }
    if (st.draft?.tool === "circle" || st.draft?.tool === "hole") {
      const o = screen(st.draft.c, w, h);
      ctx.strokeStyle = SEL;
      ctx.setLineDash([4, 3]);
      ctx.beginPath();
      ctx.arc(o.x, o.y, st.draft.r * c.k, 0, Math.PI * 2);
      ctx.stroke();
      ctx.setLineDash([]);
    }
    if (st.draft?.tool === "polygon") {
      const pts = st.draft.points;
      if (pts.length) {
        ctx.strokeStyle = SEL;
        ctx.setLineDash([4, 3]);
        ctx.beginPath();
        pts.forEach((p, i) => {
          const s = screen(p, w, h);
          if (i === 0) ctx.moveTo(s.x, s.y);
          else ctx.lineTo(s.x, s.y);
        });
        ctx.stroke();
        ctx.setLineDash([]);
        for (const p of pts) {
          const s = screen(p, w, h);
          ctx.fillStyle = SEL;
          ctx.beginPath();
          ctx.arc(s.x, s.y, 3, 0, Math.PI * 2);
          ctx.fill();
        }
      }
    }

    ctx.fillStyle = MUTED;
    ctx.font = "11px IBM Plex Mono, ui-monospace, monospace";
    ctx.fillText(`${step} mm`, 12, h - 12);

    void issues;
  }, [draft, loads, mesh, restraints, selectedNodeIds, selectedShapeId, shapes, tool, world]);

  useEffect(() => {
    let raf = 0;
    const loop = () => {
      paint();
      raf = requestAnimationFrame(loop);
    };
    raf = requestAnimationFrame(loop);
    return () => cancelAnimationFrame(raf);
  }, [paint]);

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      const t = e.target as HTMLElement | null;
      if (t && (t.tagName === "INPUT" || t.tagName === "TEXTAREA" || t.tagName === "SELECT")) return;
      const st = usePre.getState();
      if (e.key === "Escape") {
        st.setDraft(null);
        st.setTool("select");
      }
      if (e.key === "Delete" || e.key === "Backspace") st.deleteSelected();
      if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === "z") {
        e.preventDefault();
        if (e.shiftKey) st.redo();
        else st.undo();
      }
      if (e.key === "Enter" && st.draft?.tool === "polygon") {
        st.addPolygon(st.draft.points);
        st.setDraft(null);
      }
      if (e.key === "f") fit();
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [fit]);

  const onWheel = (e: React.WheelEvent) => {
    e.preventDefault();
    const rect = wrapRef.current!.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    const w = rect.width;
    const h = rect.height;
    const before = world(mx, my, w, h);
    const factor = e.deltaY < 0 ? 1.12 : 1 / 1.12;
    cam.current.k = Math.max(0.4, Math.min(80, cam.current.k * factor));
    const after = world(mx, my, w, h);
    cam.current.x += before.x - after.x;
    cam.current.y += before.y - after.y;
  };

  const pos = (e: React.PointerEvent) => {
    const rect = wrapRef.current!.getBoundingClientRect();
    return { x: e.clientX - rect.left, y: e.clientY - rect.top, w: rect.width, h: rect.height };
  };

  const onDown = (e: React.PointerEvent) => {
    (e.target as HTMLElement).setPointerCapture(e.pointerId);
    const p = pos(e);
    const wp = snapPt(world(p.x, p.y, p.w, p.h), 1);
    const st = usePre.getState();
    if (e.button === 1 || e.buttons === 4 || st.tool === "select" && e.altKey) {
      drag.current = { mode: "pan", lx: p.x, ly: p.y };
      return;
    }
    if (st.tool === "rect") {
      st.setDraft({ tool: "rect", a: wp, b: wp });
      drag.current = { mode: "draw" };
      return;
    }
    if (st.tool === "circle") {
      st.setDraft({ tool: "circle", c: wp, r: 0 });
      drag.current = { mode: "draw" };
      return;
    }
    if (st.tool === "hole") {
      st.setDraft({ tool: "hole", c: wp, r: 0 });
      drag.current = { mode: "draw" };
      return;
    }
    if (st.tool === "polygon") {
      const d = st.draft?.tool === "polygon" ? st.draft.points : [];
      if (d.length >= 3 && dist(wp, d[0]) * cam.current.k < 10) {
        st.addPolygon(d);
        st.setDraft(null);
        st.setTool("select");
        return;
      }
      st.setDraft({ tool: "polygon", points: [...d, wp] });
      return;
    }
    if (st.tool === "node" && st.mesh) {
      const n = nearestNode(st.mesh, world(p.x, p.y, p.w, p.h), 12 / cam.current.k);
      if (n) st.toggleNode(n.id, e.shiftKey);
      return;
    }
    const hit = hitShape(st.shapes, world(p.x, p.y, p.w, p.h), 6 / cam.current.k);
    st.selectShape(hit?.id ?? null);
    if (!hit && e.button === 0) drag.current = { mode: "pan", lx: p.x, ly: p.y };
  };

  const onMove = (e: React.PointerEvent) => {
    const p = pos(e);
    const wp = snapPt(world(p.x, p.y, p.w, p.h), 1);
    const d = drag.current;
    if (d?.mode === "pan") {
      cam.current.x -= (p.x - d.lx) / cam.current.k;
      cam.current.y += (p.y - d.ly) / cam.current.k;
      drag.current = { mode: "pan", lx: p.x, ly: p.y };
      return;
    }
    const st = usePre.getState();
    if (st.draft?.tool === "rect") st.setDraft({ tool: "rect", a: st.draft.a, b: wp });
    if (st.draft?.tool === "circle") st.setDraft({ tool: "circle", c: st.draft.c, r: dist(st.draft.c, wp) });
    if (st.draft?.tool === "hole") st.setDraft({ tool: "hole", c: st.draft.c, r: dist(st.draft.c, wp) });
  };

  const onUp = () => {
    const st = usePre.getState();
    const d = st.draft;
    if (d?.tool === "rect") {
      st.addRect(d.a.x, d.a.y, d.b.x - d.a.x, d.b.y - d.a.y);
      st.setDraft(null);
      st.setTool("select");
    }
    if (d?.tool === "circle") {
      st.addCircle(d.c.x, d.c.y, d.r);
      st.setDraft(null);
      st.setTool("select");
    }
    if (d?.tool === "hole") {
      st.addHole(d.c.x, d.c.y, d.r);
      st.setDraft(null);
      st.setTool("select");
    }
    drag.current = null;
  };

  return (
    <div ref={wrapRef} className="relative h-full min-h-0 w-full bg-viewport">
      <canvas
        ref={canvasRef}
        className="block h-full w-full touch-none"
        onWheel={onWheel}
        onPointerDown={onDown}
        onPointerMove={onMove}
        onPointerUp={onUp}
        onPointerCancel={onUp}
      />
      <button
        type="button"
        onClick={fit}
        className="absolute right-3 top-3 h-8 rounded-md border border-border bg-bg/90 px-2.5 text-xs text-muted hover:text-fg"
      >
        Einpassen
      </button>
    </div>
  );
}

function pathRing(
  ctx: CanvasRenderingContext2D,
  ring: Vec2[],
  w: number,
  h: number,
  screen: (p: Vec2, w: number, h: number) => { x: number; y: number },
) {
  ctx.beginPath();
  ring.forEach((p, i) => {
    const s = screen(p, w, h);
    if (i === 0) ctx.moveTo(s.x, s.y);
    else ctx.lineTo(s.x, s.y);
  });
  ctx.closePath();
}

function drawSupport(ctx: CanvasRenderingContext2D, x: number, y: number, ux: boolean, uy: boolean) {
  ctx.beginPath();
  if (uy) {
    ctx.moveTo(x, y);
    ctx.lineTo(x - 6, y + 10);
    ctx.lineTo(x + 6, y + 10);
    ctx.closePath();
    ctx.fill();
  } else if (ux) {
    ctx.moveTo(x, y);
    ctx.lineTo(x - 10, y - 6);
    ctx.lineTo(x - 10, y + 6);
    ctx.closePath();
    ctx.fill();
  }
}

function drawArrow(ctx: CanvasRenderingContext2D, x: number, y: number, dx: number, dy: number) {
  ctx.beginPath();
  ctx.moveTo(x, y);
  ctx.lineTo(x + dx, y + dy);
  ctx.stroke();
  const a = Math.atan2(dy, dx);
  ctx.beginPath();
  ctx.moveTo(x + dx, y + dy);
  ctx.lineTo(x + dx - 7 * Math.cos(a - 0.4), y + dy - 7 * Math.sin(a - 0.4));
  ctx.lineTo(x + dx - 7 * Math.cos(a + 0.4), y + dy - 7 * Math.sin(a + 0.4));
  ctx.closePath();
  ctx.fill();
}

function dimLabel(
  ctx: CanvasRenderingContext2D,
  a: Vec2,
  b: Vec2,
  w: number,
  h: number,
  screen: (p: Vec2, w: number, h: number) => { x: number; y: number },
) {
  const mid = screen({ x: (a.x + b.x) / 2, y: (a.y + b.y) / 2 }, w, h);
  ctx.fillStyle = INK;
  ctx.font = "11px IBM Plex Mono, ui-monospace, monospace";
  ctx.fillText(`${Math.abs(b.x - a.x).toFixed(0)} × ${Math.abs(b.y - a.y).toFixed(0)} mm`, mid.x + 8, mid.y - 8);
  void DANGER;
}
