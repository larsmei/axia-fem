declare module "@/wasm/axia_fem.js" {
  export function preview_inp(inp: string): string;
  export function solve_inp(inp: string): string;
  export default function init(opts?: {
    module_or_path?: string | URL | Request | BufferSource;
  }): Promise<unknown>;
}

declare module "*.wasm?url" {
  const src: string;
  export default src;
}
