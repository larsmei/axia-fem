# Axia

Linear-statischer FEM-Solver im Browser. CalculiX-kompatibles INP einlesen, rechnen, verformtes Netz und Spannungsfelder anzeigen, FRD exportieren.

Der Solver ist in Rust geschrieben, als WebAssembly gebaut und läuft komplett clientseitig.

## Features

- Elemente: **C3D8**, **C3D4**, **CPS4**, **CPE4**, **CPS3**, **CPE3**
- INP-Parser (Knoten, Elemente, Material, Randbedingungen, Lasten)
- Linear-statische Analyse (isotroper Elastizität)
- 3D-Viewer: undeformiert / deformiert, Felder (von Mises, |u|, ux/uy/uz, Sij)
- FRD-Export
- Mitgelieferte Beispiele (Balken, Platte, Hex-Gitter)

## Entwicklung

Voraussetzungen: Node.js 22+, npm.

```bash
npm install
npm run dev
```

Die App läuft dann lokal (Port 8080).

```bash
npm run build
npm run typecheck
```

## WASM neu bauen

Rust-Toolchain mit `wasm32-unknown-unknown` und `wasm-bindgen-cli` (Version 0.2.128):

```bash
cd solver
cargo build --release --target wasm32-unknown-unknown
wasm-bindgen --target web --out-dir ../src/wasm \
  target/wasm32-unknown-unknown/release/axia_fem.wasm
cp ../src/wasm/axia_fem_bg.wasm ../public/axia_fem_bg.wasm
```

## Lizenz

MIT — siehe [LICENSE](LICENSE).
