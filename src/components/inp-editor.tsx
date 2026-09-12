import { useEffect, useMemo, useRef } from "react";
import { cn } from "@/lib/utils";

type Props = {
  value: string;
  onChange: (v: string) => void;
  className?: string;
};

export function InpEditor({ value, onChange, className }: Props) {
  const preRef = useRef<HTMLPreElement>(null);
  const taRef = useRef<HTMLTextAreaElement>(null);
  const lines = useMemo(() => value.split("\n"), [value]);

  useEffect(() => {
    const ta = taRef.current;
    const pre = preRef.current;
    if (!ta || !pre) return;
    const sync = () => {
      pre.scrollTop = ta.scrollTop;
      pre.scrollLeft = ta.scrollLeft;
    };
    ta.addEventListener("scroll", sync);
    return () => ta.removeEventListener("scroll", sync);
  }, []);

  return (
    <div className={cn("relative min-h-0 flex-1 overflow-hidden bg-code", className)}>
      <pre
        ref={preRef}
        className="pointer-events-none absolute inset-0 overflow-hidden font-mono text-[12.5px] leading-5 text-fg"
      >
        <code className="block py-3">
          {lines.map((line, i) => (
            <span key={i} className="flex">
              <span className="w-11 shrink-0 pr-2 text-right text-[11px] text-subtle">
                {i + 1}
              </span>
              <span className="min-w-0 flex-1 pr-3 whitespace-pre">
                {highlight(line)}
              </span>
            </span>
          ))}
        </code>
      </pre>
      <textarea
        ref={taRef}
        value={value}
        spellCheck={false}
        onChange={(e) => onChange(e.target.value)}
        className="absolute inset-0 resize-none overflow-auto bg-transparent py-3 pr-3 pl-11 font-mono text-[12.5px] leading-5 text-transparent caret-fg outline-none"
        aria-label="CalculiX INP"
      />
    </div>
  );
}

function highlight(line: string) {
  const t = line.trimStart();
  if (t.startsWith("**")) {
    return <span className="text-subtle">{line || " "}</span>;
  }
  if (t.startsWith("*")) {
    const idx = line.indexOf("*");
    const rest = line.slice(idx);
    const comma = rest.indexOf(",");
    const kw = comma === -1 ? rest : rest.slice(0, comma);
    const params = comma === -1 ? "" : rest.slice(comma);
    return (
      <>
        {line.slice(0, idx)}
        <span className="text-keyword">{kw}</span>
        <span className="text-param">{params}</span>
      </>
    );
  }
  return <span className="text-muted">{line || " "}</span>;
}
