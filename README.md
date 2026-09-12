# Axia

Linear-statischer FEM-Solver. CalculiX-/Abaqus-kompatibles INP einlesen, rechnen, FRD/DAT schreiben.

Zwei Frontends, ein Solver:

- **CLI** (`axia`) — natives Binary für Linux, Windows und macOS
- **Browser** — Rust als WebAssembly, 3D-Viewer, komplett clientseitig

**Repo:** [larsmei/axia-fem](https://github.com/larsmei/axia-fem) · **Releases:** [latest](https://github.com/larsmei/axia-fem/releases)

## CLI

Wie CalculiX: Job-Name ohne Endung. `axia job` liest `job.inp` und schreibt `job.frd` plus `job.dat`.

```bash
axia --version
axia --help

axia examples/patch_c3d8          # → examples/patch_c3d8.frd / .dat
axia examples/cantilever_b32.inp
axia --check model.inp            # nur parsen
axia --json job.inp               # Statistik als JSON
cat deck.inp | axia - --stdout > out.frd
```

Mitgelieferte Decks in [`examples/`](examples/):

| Datei | Inhalt |
|---|---|
| `patch_c3d8.inp` | Zug-Patch C3D8, \(u_x = 0{,}01\) |
| `cantilever_b32.inp` | Timoshenko-Kragträger B32 |
| `plate_s4r.inp` | gelenkig gelagerte S4R-Platte |

### Binary bauen

```bash
cd solver
cargo build --release --bin axia
./target/release/axia --version
```

Cross-Compile (Linux-Host) und Packen:

```bash
# Targets + Linker (Debian/Ubuntu)
sudo apt-get install -y gcc-mingw-w64-x86-64 gcc-aarch64-linux-gnu musl-tools
rustup target add x86_64-unknown-linux-musl aarch64-unknown-linux-gnu x86_64-pc-windows-gnu

./solver/scripts/release.sh
# → dist/axia-<ver>-<triple>.tar.gz  /  .zip
```

GitHub Actions (`.github/workflows/release.yml`) baut bei einem Tag `v*` zusätzlich **macOS** (Intel + Apple Silicon) und **Windows MSVC**.

## Features

- INP-Parser: `*NODE`, `*ELEMENT`, `*NSET`/`*ELSET` (+ `GENERATE`), `*MATERIAL`/`*ELASTIC`/`*DENSITY`
- Schnitte: `*SOLID SECTION`, `*SHELL SECTION`, `*BEAM SECTION` (`RECT`, `CIRC`, `PIPE`, `GENERAL`)
- Lasten und Lager: `*BOUNDARY` (DOF 1–6), `*CLOAD`, `*DLOAD` (`P`, `P1…P6`, `GRAV`, `PX`/`PY`/`PZ`)
- Linear-statische Analyse, isotrope Elastizität
- Sparse-Assembly (Triplet → CSR), Cholesky oder PCG
- FRD- und DAT-Export (cgx-kompatible Elementtypen)
- 3D-Viewer (Three.js): undeformiert / deformiert, von Mises, |u|, ux/uy/uz, Sij

## Elementbibliothek

### Kontinuum

| Typ | Knoten | DOF/Knoten | Bemerkung |
|---|---|---|---|
| C3D8 / C3D8R | 8 | 3 | Hexaeder, 2×2×2 |
| C3D20 / C3D20R | 20 | 3 | quadratisches Serendipity-Hexaeder, 3×3×3 / 2×2×2 |
| C3D4 | 4 | 3 | Tetraeder |
| C3D10 | 10 | 3 | quadratisches Tetraeder, 4-Punkt |
| CPS4 / CPE4 | 4 | 2 | Scheibe, Spannungs-/Dehnungszustand |
| CPS8 / CPE8 / CPS8R | 8 | 2 | quadratische Scheibe, 3×3 / 2×2 |
| CPS3 / CPE3 | 3 | 2 | Dreiecksscheibe |
| CPS6 / CPE6 | 6 | 2 | quadratisches Dreieck, 3-Punkt |

Knotenreihenfolge der quadratischen Elemente wie Abaqus/CalculiX: **Ecken zuerst, dann Kantenmittelpunkte**. Serendipity-Netze haben keine Flächenmittelpunkte.

### Schale (Reissner–Mindlin, 6 DOF)

| Typ | Knoten | Formulierung |
|---|---|---|
| S4 / S4R | 4 | MITC4 — kein Schub-Locking |
| S8 / S8R | 8 | quadratisch, selektiv reduzierter Schub |
| S3 / S3R | 3 | DKT-Biegung + CST-Membran |
| S6 | 6 | quadratisches Dreieck |

Dicke über `*SHELL SECTION`. Positive Drucklast `P` wirkt in Richtung der positiven Elementnormalen (rechte Hand, Knoten CCW).

```
*ELEMENT, TYPE=S4R, ELSET=PLATE
1, 1, 2, 3, 4
*SHELL SECTION, ELSET=PLATE, MATERIAL=STEEL
1.0
*BOUNDARY
1, 1, 6
*DLOAD
EALL, P, -0.01
```

### Balken (Timoshenko, 6 DOF)

| Typ | Knoten | Reihenfolge |
|---|---|---|
| B31 | 2 | Ende 1, Ende 2 |
| B32 | 3 | **Ende 1, Ende 2, Mitte** |

Querschnitt: `RECT`, `CIRC`, `PIPE`, `GENERAL`. Die zweite Datenzeile von `*BEAM SECTION` ist die n1-Richtung.

```
*ELEMENT, TYPE=B32, ELSET=BEAM
1, 1, 3, 2
*BEAM SECTION, ELSET=BEAM, MATERIAL=STEEL, SECTION=RECT
10, 20
0, 0, 1
*BOUNDARY
1, 1, 6
*CLOAD
3, 2, -100
```

Einspannung eines Kragträgers: `1, 1, 6` (Translation + Rotation). `PX`/`PY`/`PZ` sind Streckenlasten.

Gemischte Modelle (Kontinuum + Schale/Balken) verwenden 6 DOF pro Knoten. Unbenutzte Rotationen an reinen Kontinuumsknoten werden automatisch festgehalten.

## Beispiele und Checks

| Beispiel | Erwartung |
|---|---|
| Zugstab C3D8 / C3D20 | Patch-Test: \(u_x = 0{,}01\), \(\sigma_{xx} = 210\) |
| Kragträger B32 | Timoshenko \(\delta = PL^3/3EI + PL/kAG \approx 0{,}1906\) |
| Kragplatte S4R | Euler \(\delta \approx 1{,}905\) (8×2, \(\nu=0\): \(\lvert u\rvert_{\max} \approx 1{,}83\)) |
| Quadratplatte S4R | Kirchhoff \(\delta_{\max} \approx 0{,}00406\,qa^4/D \approx 0{,}211\) |
| Kragträger CPS8 | Euler \(\delta \approx 1{,}905\), näher als lineare CPS4 |

```bash
cd solver && cargo test
./target/release/axia ../examples/patch_c3d8.inp --json
```

## Web-App

Voraussetzungen: Node.js 22+, npm, Rust mit Target `wasm32-unknown-unknown`, `wasm-bindgen-cli` 0.2.128.

```bash
npm install
npm run dev
```

```bash
npm run build
npm run typecheck
```

### WASM neu bauen

```bash
cd solver
cargo build --release --target wasm32-unknown-unknown --lib --no-default-features
wasm-bindgen --target web --out-dir ../src/wasm \
  target/wasm32-unknown-unknown/release/axia_fem.wasm
cp ../src/wasm/axia_fem_bg.wasm ../public/axia_fem_bg.wasm
```

## Lizenz

MIT — siehe [LICENSE](LICENSE).
