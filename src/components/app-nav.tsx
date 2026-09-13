import { Link } from "@tanstack/react-router";
import { cn } from "@/lib/utils";

export function AppMark() {
  return (
    <svg width="28" height="28" viewBox="0 0 28 28" aria-hidden className="shrink-0">
      <rect width="28" height="28" rx="7" className="fill-surface-2 stroke-border" strokeWidth="1" />
      <path d="M8 20 L14 7 L20 20 Z" fill="none" className="stroke-fg" strokeWidth="1.4" strokeLinejoin="round" />
      <circle cx="14" cy="7" r="1.4" className="fill-fg" />
      <circle cx="8" cy="20" r="1.4" className="fill-fg" />
      <circle cx="20" cy="20" r="1.4" className="fill-fg" />
    </svg>
  );
}

export function AppNav({ active }: { active: "solver" | "pre" }) {
  const item = (to: "/" | "/pre", label: string, on: boolean) => (
    <Link
      to={to}
      className={cn(
        "inline-flex h-8 items-center rounded-md px-2.5 text-xs font-medium",
        on ? "bg-surface text-fg" : "text-muted hover:text-fg",
      )}
    >
      {label}
    </Link>
  );
  return (
    <nav className="flex rounded-md border border-border bg-surface-2 p-0.5" aria-label="Bereich">
      {item("/", "Solver", active === "solver")}
      {item("/pre", "Präprozessor", active === "pre")}
    </nav>
  );
}
