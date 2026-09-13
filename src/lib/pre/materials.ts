import type { Material } from "./types";

export const CATALOG: Material[] = [
  { id: "steel", name: "Stahl S235", E: 210000, nu: 0.3, density: 7.85e-9, thickness: 1 },
  { id: "alu", name: "Aluminium 6061", E: 70000, nu: 0.33, density: 2.7e-9, thickness: 1 },
  { id: "concrete", name: "Beton C30", E: 32000, nu: 0.2, density: 2.4e-9, thickness: 10 },
];

export function defaultMaterials(): Material[] {
  return CATALOG.map((m) => ({ ...m }));
}
