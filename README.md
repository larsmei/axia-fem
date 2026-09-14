# Axia

FEM-Solver. CalculiX-/Abaqus-kompatibles INP einlesen, rechnen, FRD/DAT schreiben.

Zwei Frontends, ein Solver:

- **CLI** (`axia`) — natives Binary für Linux, Windows und macOS
- **Browser** — Rust als WebAssembly, 3D-Viewer, komplett clientseitig
- **Präprozessor** (ab 1.11) — 3D-Volumen (Quader/Zylinder/Kugel, C3D8/C3D6) und 2D-Scheibe, Vernetzung, Material, Flächenlager/Lasten, INP-Export

**Repo:** [larsmei/axia-fem](https://github.com/larsmei/axia-fem) · **Releases:** [v1.11.0](https://github.com/larsmei/axia-fem/releases/tag/v1.11.0)

## Downloads (CLI 1.11.0)

Fertige Binaries am Release. Archiv enthält `axia` / `axia.exe`, `README`, `LICENSE` und `examples/`.

| Plattform | Datei |
|---|---|
| Windows x64 | [axia-1.11.0-windows-x64.zip](https://github.com/larsmei/axia-fem/releases/download/v1.11.0/axia-1.11.0-windows-x64.zip) |
| Linux x64 (glibc) | [axia-1.11.0-linux-x64.tar.gz](https://github.com/larsmei/axia-fem/releases/download/v1.11.0/axia-1.11.0-linux-x64.tar.gz) |
| Linux x64 (musl, statisch) | [axia-1.11.0-linux-x64-musl.tar.gz](https://github.com/larsmei/axia-fem/releases/download/v1.11.0/axia-1.11.0-linux-x64-musl.tar.gz) |
| Linux ARM64 | [axia-1.11.0-linux-arm64.tar.gz](https://github.com/larsmei/axia-fem/releases/download/v1.11.0/axia-1.11.0-linux-arm64.tar.gz) |
| macOS Apple Silicon | [axia-1.11.0-macos-arm64.tar.gz](https://github.com/larsmei/axia-fem/releases/download/v1.11.0/axia-1.11.0-macos-arm64.tar.gz) |
| macOS Intel | [axia-1.11.0-macos-x64.tar.gz](https://github.com/larsmei/axia-fem/releases/download/v1.11.0/axia-1.11.0-macos-x64.tar.gz) |

Prüfsummen: [SHA256SUMS](https://github.com/larsmei/axia-fem/releases/download/v1.11.0/SHA256SUMS). Standard-Solver ist **faer** (kein Intel-MKL nötig).

**1.11.0** — Präprozessor 3D (Quader/Zylinder/Kugel, C3D8/C3D6, Flächenlager, Orbit-Viewport) und breitere CalculiX-INP-Kompatibilität (Fortran-`1.d0`, verschachtelte Sets, `CENTRIF`, `TYPE=ORTHO`/`ENGINEERING`, `BEAM SECTION=BOX`, `TRANSFORM TYPE=C`, Composite-Schalen, `ROT NODE`, lineare Fallbacks für NLGEOM/Kontakt). Vergleichsläufer `solver/scripts/ccx-compare.py` gegen die ccx-2.23-Suite. **1.10.0** — Browser-Präprozessor: Rechteck/Kreis/Polygon, strukturiertes Quad-Netz und Delaunay, Stahl/Alu/Beton, Kantenlager und Kräfte, Modellcheck, CalculiX-INP (CPS3/CPS4), Übergabe an den WASM-Solver. Vorlagen Kragarm, gelochte Platte, Zugstab. **1.9.8** — `*TIE`/`*COUPLING`/`*CONTACT` auf C3D20/C3D10/C3D15-Flächen inkl. Mittelknoten (vorher nur Ecken → klaffende Kanten). **1.9.7** — MKL: kein stiller Prozessabbruch mehr. Unvollständige Redist (ohne `mkl_core` / CPU-Kernel `mkl_avx2`/`mkl_def`) wird vor `pardiso()` erkannt; Selbsttest in einem Kindprozess fängt Intel-`abort()` (FATAL ERROR / OpenMP #15). `--solver auto` fällt dann auf **faer** zurück, `--solver mkl` schreibt die Ursache auf stderr. `KMP_DUPLICATE_LIB_OK=TRUE`, `perm[n]`, Matrix-Checker. **1.9.6** — FRD: C3D20/C3D20R und C3D15 mit ccx-Midside-Permutation (Abaqus-INP → FAM/he20), sonst verdrehtes Hex20 in Mecway/cgx. **1.9.5** — FRD: `MASS`/`ROTARYI` wie ccx nicht ins Elementnetz (Mecway-Typ 11 braucht zwei I10-Felder; sonst „Field of length 10 missing at column 14“). **1.9.4** — Sparse-Solver per `--solver` / `AXIA_SOLVER`: **faer** (reines Rust, supernodales \(LL^\top\)/\(LU\), keine extra Library, PARDISO-Klasse für SPD). `auto` = MKL → Panua → faer → rivrs-sparse. **1.9.3** — Mecway/CalculiX: `*ELEMENT` ohne `ELSET=` bekommt Material/Dicke/`*MASS` über spätere `*ELSET`-Mitgliedschaft; MASS/DASHPOT/GAP brauchen kein Kontinuum-`*MATERIAL`; GRAV auf `*MASS` als \(F=mg\). **1.9.2** — Windows-MKL: `mkl_rt.dll` **neben `axia.exe`** (Shim `libmkl_rt.dll` im selben Ordner, nicht in `%TEMP%`). **1.9.1** — FRD long-ASCII wie ccx (Spalte 74, Fortran `E-02`, ` -1`/` -2`/` -3`), Windows-MKL findet `mkl_rt.dll` unter `$MKLROOT/bin`, Viewer Jet-Farbverlauf / Mecway-Colorbar. **1.9** — fehlende Elemente: CAX3/CAX6, C3D10T (B-bar), MASS, ROTARYI, DASHPOTA, GAPUNI, Aliase B21/B22 und C3D20RI. **1.8** — `*STATIC, RIKS` (Crisfield-Bogenlänge, Snap-Through). **1.7.1** — Viewer zeichnet T3D2/T3D3/SPRINGA (Fachwerk NLGEOM). **1.7** — Coulomb-`*FRICTION`, `NLGEOM`+`*PLASTIC` auf Kontinuum. **1.6** — `*CONTACT PAIR` Node-to-Surface, Penalty, reibungsfrei. **1.5** — `*PLASTIC` J2 für Kontinuum (C3D*), PEEQ im FRD. **1.4** — `*STEP, NLGEOM` für Kontinuum (C3D8/20/4/10/6/15), Total-Lagrange St. Venant–Kirchhoff. **1.3** — mehrere `*STEP`, `*CONTROLS`. **1.2** — C3D15, CAX, Membran, Kontinuum-Beulen, Wärme+. **1.1** — Wärme, Dynamik, NLGEOM/`*PLASTIC` (T3D2).

## Präprozessor (Browser)

Unter **Präprozessor** (Route `/pre`) entsteht ein 3D-Volumenmodell oder ein 2D-Scheibenmodell, das der Solver direkt rechnet.

1. Geometrie — 3D: Quader, Zylinder, Kugel, extrudiertes Polygon; 2D: Rechteck, Kreis, Polygon, Löcher
2. Netz — C3D8-Hexeder / C3D6-Prismen (3D) oder CPS4/CPS3 (2D), Kantenlänge, Qualitätscheck
3. Material — Stahl, Aluminium, Beton oder eigene Kennwerte (E, ν, ρ; Dicke nur 2D)
4. Randbedingungen — Lager und Kräfte an Flächen (X± Y± Z±) oder Knoten, optional Eigengewicht
5. Prüfen — fehlendes Netz/Material/Lager, Mechanismus
6. Export — CalculiX-INP herunterladen oder **Lösen** an den WASM-Solver übergeben

Einheiten: **mm, N, MPa** (Dichte in t/mm³). Vorlagen: Kragarm, Zugstab, Druckwürfel, Zylinder, gelochte Platte.

Viewport: Orbit (links ziehen), Schieben (mittlere Taste), Zoom (Rad). Iso / Draufsicht / Einpassen. Z zeigt nach oben.

## CLI

Wie CalculiX: Job-Name ohne Endung. `axia job` liest `job.inp` und schreibt `job.frd` plus `job.dat`.

```bash
axia --version
axia --help

axia examples/patch_c3d8          # → examples/patch_c3d8.frd / .dat
axia examples/cantilever_b32.inp
axia --check model.inp            # nur parsen
axia --json job.inp               # Statistik als JSON
axia --solver faer job.inp        # reines Rust, keine extra Library
cat deck.inp | axia - --stdout > out.frd
```

Mitgelieferte Decks in [`examples/`](examples/):

| Datei | Inhalt |
|---|---|
| `patch_c3d8.inp` | Zug-Patch C3D8, \(u_x = 0{,}01\) |
| `patch_c3d6.inp` | Zug-Patch C3D6-Wedge |
| `patch_c3d15.inp` | C3D15, confined \(\varepsilon_z=0{,}001\) |
| `pipe_cax4.inp` | dickwandiges Rohr CAX4, Innendruck (Lamé) |
| `pipe_cax3.inp` | dasselbe mit CAX3-Dreiecken |
| `patch_c3d10t.inp` | C3D10T B-bar, confined \(\varepsilon_x=0{,}001\) |
| `mass_sdof.inp` | `*MASS` + SPRINGA, \(f\approx 1{,}59\,\mathrm{Hz}\) |
| `gapuni.inp` | GAPUNI geschlossen, \(u=0{,}01\) |
| `dashpota.inp` | DASHPOTA überdämpft, \(\lvert u(T)\rvert\ll u_0\) |
| `membrane_m3d4.inp` | M3D4-Membran-Patch, \(u_x=0{,}01\) |
| `cantilever_b32.inp` | Timoshenko-Kragträger B32 |
| `plate_s4r.inp` | gelenkig gelagerte S4R-Platte |
| `truss_t3d2.inp` | Fachwerkstab T3D2, \(u = FL/EA\) |
| `equation_bars.inp` | zwei Stäbe gekoppelt mit `*EQUATION` |
| `tie_bars.inp` | dasselbe mit `*TIE` |
| `freq_t3d2.inp` | axiale Eigenfrequenz T3D2 |
| `buckle_b31.inp` | Euler-Beulen Kragträger B31 |
| `include_main.inp` | `*INCLUDE` (zieht `include_mat.inp`) |
| `thermal_bar.inp` | eingespannter Stab, `*EXPANSION` / `*TEMPERATURE` |
| `nlgeom_truss.inp` | T3D2 `*STEP, NLGEOM`, kleine Dehnung \(u=FL/EA\) |
| `nlgeom_c3d8.inp` | C3D8 `*STEP, NLGEOM`, Patch \(u_x=0{,}01\) |
| `nlgeom_stretch.inp` | C3D8 SVK, \(\lambda=1{,}2\), Cauchy \(\sigma_{xx}=26{,}4\) |
| `plastic_bar.inp` | T3D2 `*PLASTIC` (isotrope Verfestigung), \(u\approx 3{,}095\) |
| `plastic_c3d8.inp` | C3D8 J2, \(\sigma=250\), \(u_x\approx 0{,}03095\) |
| `contact_blocks.inp` | zwei C3D8, `*CONTACT PAIR`, \(u_{z,\mathrm{mid}}=-0{,}005\) |
| `contact_friction.inp` | Coulomb-Haften \(\mu=0{,}8\), \(u_{x,\mathrm{slave}}\approx 0\) |
| `nlgeom_plastic_c3d8.inp` | C3D8 `NLGEOM`+`*PLASTIC`, \(u_x\approx 0{,}03095\) |
| `riks_truss.inp` | von Mises-Fachwerk, `*STATIC, RIKS`, Snap-Through |
| `riks_c3d8.inp` | C3D8 Riks, Patch \(u_x=0{,}01\), \(\lambda=1\) |
| `heat_bar.inp` | stationäre Wärmeleitung T3D2, \(T(L/2)=50\) |
| `dynamic_sdof.inp` | Newmark, SDOF \(u(T/2)=-u_0\) |

### Binary bauen

```bash
cd solver
cargo build --release --bin axia
./target/release/axia --version
```

### Sparse-Solver (nativ)

Auswahl per CLI `--solver` / `-s` oder Umgebungsvariable `AXIA_SOLVER`:

| Name | Was | Extra-Library |
|---|---|---|
| `auto` (Standard) | MKL-PARDISO → Panua-PARDISO → **faer** → rivrs-sparse → dens/PCG | nur für PARDISO |
| `mkl` | Intel MKL PARDISO | `mkl_rt` (+ Kern-DLLs) |
| `panua` | Panua PARDISO | `libpardiso` |
| `pardiso` | MKL oder Panua, ohne Rust-Fallback | wie oben |
| `faer` | supernodales \(LL^\top\), sonst \(LU\) | **keine** — reines Rust |
| `rivrs` | rivrs-sparse \(LDL^\top\) (APTP, METIS/AMD) | **keine** (METIS ist optional einkompiliert) |
| `cholesky` | dichte In-Crate-Cholesky | — (\(n\le 900\)) |
| `pcg` | vorkonditioniertes CG | — |

`faer` ist der nächste reine-Rust-Solver an PARDISO: supernodale Cholesky-Faktorisierung (wie CHOLMOD), AMD-Ordering, parallele BLAS-3-Kerne, bereits Abhängigkeit von Axia. Für SPD-Steifigkeitsmatrizen (linear-statisch nach Randbedingungen) ist das der richtige Algorithmus. Bei indefiniter \(K\) (selten) fällt faer intern auf supernodales \(LU\) zurück.

```bash
axia --solver faer job.inp
axia --solver rivrs job.inp
AXIA_SOLVER=faer axia job
```

Beim Start steht auf stderr z. B. `axia: sparse solver: faer (supernodal LLT)`. Die JSON-Statistik enthält dasselbe Feld `solver`.

Die Browser-WASM-Variante verwendet weiterhin die eingebaute Cholesky-/PCG-Kette (kein PARDISO / faer im Browser).

#### Intel MKL unter Windows

`pardiso-wrapper` 0.1.2 sucht nach **`libmkl_rt.dll`**. Intel oneAPI liefert **`mkl_rt.dll`**. Axia sucht in dieser Reihenfolge:

1. **Ordner von `axia.exe`** (portable: DLLs neben die EXE legen)
2. aktuelles Arbeitsverzeichnis
3. `$MKLROOT\bin`, `bin\intel64`, `redist\intel64`, `lib`, …
4. `PATH`

und legt bei Bedarf **`libmkl_rt.dll` im selben Ordner** wie `mkl_rt.dll` an (Hardlink/Kopie). Der Shim darf nicht in `%TEMP%` allein liegen: `mkl_rt` lädt `mkl_core` aus **seinem eigenen Verzeichnis**.

**Nur `mkl_rt.dll` + `libiomp5md.dll` reicht nicht.** Mindestens:

| Datei | Rolle |
|---|---|
| `mkl_rt.dll` / `mkl_rt.2.dll` | Dispatcher (SDL) |
| `libmkl_rt.dll` | Name, den der Wrapper sucht (Axia erzeugt ihn) |
| `mkl_core.2.dll` | Kern |
| `mkl_intel_thread.2.dll` oder `mkl_sequential.2.dll` | Threading |
| `mkl_avx2.2.dll` oder `mkl_def.2.dll` | CPU-Kernels |
| `libiomp5md.dll` | Intel OpenMP (bei `intel_thread`) |

Einfach den Inhalt von `%MKLROOT%\bin` (plus `libiomp5md.dll` aus `compiler\latest\bin`) **neben `axia.exe` kopieren**, oder:

```bat
set MKLROOT=C:\Program Files (x86)\Intel\oneAPI\mkl\latest
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
axia job
```

Wenn PARDISO nicht lädt **oder beim Selbsttest abstürzt**, schreibt Axia die Ursache auf stderr (fehlende DLLs, `Intel MKL FATAL ERROR`, OpenMP Error #15). `--solver auto` fällt danach auf **faer** zurück statt den Prozess zu beenden. `--solver mkl` bricht mit dieser Meldung ab, nicht still.

Nur `mkl_rt.dll` reicht nicht. Fehlen CPU-Kernel (`mkl_avx2.2.dll` / `mkl_def.2.dll`), ruft Intel intern `abort()` auf — bis 1.9.6 ohne Axia-Text, in Mecway unsichtbar. Ab 1.9.7 prüft Axia die Dateien **bevor** PARDISO läuft und fängt einen verbleibenden C-`abort()` in einem Kindprozess.

Umgebungsvariablen (optional): `KMP_DUPLICATE_LIB_OK=TRUE` (gesetzt, wenn nicht vorhanden), `MKL_THREADING_LAYER=SEQUENTIAL` oder `INTEL`, `AXIA_SKIP_MKL_PROBE=1` (Selbsttest überspringen).

#### FRD (CalculiX / Mecway / cgx)

Ausgabe ist long-ASCII wie ccx `frd.c` / `frdheader.c`:

- 2C/3C-Zeilen genau **74 Zeichen**, Formatflag `1` in Spalte 74
- Datensätze ` -1` / ` -2` / ` -3` (führendes Leerzeichen)
- C3D20 / C3D15: Midside-Knoten wie ccx `frd.c` (FAM/he20, nicht Abaqus-INP)
- MASS / ROTARYI nicht im 3C-Elementblock (wie ccx; FRD-Typ 11 braucht zwei Knoten)
- Fortran **ES12.5** (` 6.00000E-02`, nicht Rust `E-2`)
- 100CL 75 Zeichen mit Flag in Spalte 75; 1PSTEP 70 Zeichen

Damit lesbar in Mecway, cgx und FreeCAD.

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

- INP-Parser: `*NODE`, `*ELEMENT`, `*NSET`/`*ELSET` (+ `GENERATE`), `*MATERIAL`/`*ELASTIC`/`*DENSITY`, `*INCLUDE`, `*EQUATION`, `*SURFACE`
- Schnitte über ELSET-Mitgliedschaft (CalculiX/Mecway): `*ELEMENT` darf `ELSET=` weglassen; `*SHELL SECTION` / `*SOLID SECTION` / `*MASS` gelten für jedes Element, das später in dem Set steht
- Schnitte: `*SOLID SECTION`, `*SHELL SECTION`, `*MEMBRANE SECTION`, `*BEAM SECTION` (`RECT`, `CIRC`, `PIPE`, `GENERAL`), `*SPRING`, `*MASS`, `*ROTARY INERTIA`, `*DASHPOT`, `*GAP`
- Lasten und Lager: `*BOUNDARY` (DOF 1–6, NT=11), `*CLOAD`, `*DLOAD` (`P`, `P1…P6`, `GRAV`, `PX`/`PY`/`PZ`)
- Linear-statische Analyse, isotrope Elastizität, MPC-Elimination (`*EQUATION`)
- `*FREQUENCY` (lumped mass, inverse subspace) und `*BUCKLE` (geometrische Steifigkeit für T3D2, B31/B32 **und** C3D8/20/10/4/6/15)
- `*HEAT TRANSFER, STEADY STATE` — T3D2, B31, C3D8, C3D20, C3D10, C3D4, C3D6, C3D15, CPS4, CPS8, S4; `*CONDUCTIVITY`, `*DFLUX`/`*FILM`/`*CFLUX`, NT = DOF 11
- `*DYNAMIC` — implizites Newmark (\(\beta=1/4,\gamma=1/2\)), `*DAMPING` (Rayleigh), `*AMPLITUDE`, `*INITIAL CONDITIONS`
- `*EXPANSION` + `*TEMPERATURE` (isotrope Wärmedehnung, T3D2 und C3D8*)
- `*STEP, NLGEOM` — geometrisch nichtlineare Statik
  - **T3D2** korotational, Newton (`*CONTROLS, MAXITER=`, `RTOL=`)
  - **Kontinuum** C3D8/C3D8I/C3D8R/C3D20/C3D20R/C3D4/C3D10/C3D6/C3D15: Total-Lagrange, St. Venant–Kirchhoff. Kleine Dehnung reproduziert den linearen Patch-Test. Lasten auf der Referenzkonfiguration (dead load).
- `*STATIC, RIKS` — modifiziertes Riks-/Crisfield-Bogenlängenverfahren (Snap-Through, Lastfaktor \(\lambda\)). Schaltet NLGEOM ein. T3D2 und Kontinuum.
- Mehrere `*STEP` / `*END STEP` — Lasten und Lager kumulativ, letzter Schritt bestimmt die Ausgabe
- `*PLASTIC` — J2 mit isotroper Verfestigung (Kurve \(\sigma_y(\bar\varepsilon^p)\))
  - **T3D2** 1D-Return-Map
  - **Kontinuum** C3D8/C3D8I/C3D8R/C3D20/C3D20R/C3D4/C3D10/C3D6/C3D15: Radial-Return, konsistente Tangente, FRD-Block `PEEQ`
  - `NLGEOM` + `*PLASTIC` auf Kontinuum: J2 auf Green–Lagrange / PK2
- `*CONTACT PAIR` — Penalty Node-to-Surface, `*FRICTION` Coulomb (small sliding)
- Sparse-Assembly (Triplet → CSR)
- Native Sparse-Solver: **PARDISO** (Intel MKL / Panua), **faer** (reines Rust, supernodal), **rivrs-sparse**; wählbar mit `--solver` / `AXIA_SOLVER`
- Der jeweils verwendete Solver wird beim Aufruf ausgegeben (`axia: sparse solver: …`) und steht in der Statistik
- FRD- und DAT-Export (ccx/Mecway/cgx long-ASCII, Spalte 74)
- 3D-Viewer (Three.js): undeformiert / deformiert, Jet-Farbverlauf (von Mises, |u|, ux/uy/uz, Sij), Mecway-Colorbar, Achsenkreuz

## Elementbibliothek

### Kontinuum

| Typ | Knoten | DOF/Knoten | Bemerkung |
|---|---|---|---|
| C3D8 | 8 | 3 | Hexaeder, volle 2×2×2 Integration |
| C3D8I | 8 | 3 | inkompatible Wilson/Taylor-Moden, kein Schub-Locking |
| C3D8R | 8 | 3 | 1-Punkt + Hourglass-Stabilisierung |
| C3D20 / C3D20R / C3D20RI | 20 | 3 | quadratisches Serendipity-Hexaeder, 3×3×3 / 2×2×2 |
| C3D6 | 6 | 3 | linearer Wedge/Pentaeder |
| C3D15 | 15 | 3 | quadratischer Wedge, 3×3 Gauss (Dreieck × ζ) |
| C3D4 | 4 | 3 | Tetraeder |
| C3D10 | 10 | 3 | quadratisches Tetraeder, 4-Punkt |
| C3D10T | 10 | 3 | B-bar / mittlere Dilatation (weniger Volumen-Locking) |
| CAX4 / CAX4R | 4 | 2 | Achsensymmetrie \(r,z\); Gewicht \(2\pi r\) |
| CAX8 / CAX8R | 8 | 2 | quadratisch achsensymmetrisch |
| CAX3 | 3 | 2 | lineares Achsensymmetrie-Dreieck |
| CAX6 | 6 | 2 | quadratisches Achsensymmetrie-Dreieck |
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

Dicke über `*SHELL SECTION`. Positive Drucklast `P` wirkt in Richtung der positiven Elementnormalen (rechte Hand, Knoten CCW). Mecway schreibt oft `*ELEMENT, TYPE=S4` **ohne** `ELSET=` und weist das Material erst über ein späteres `*ELSET` plus `*SHELL SECTION` zu — Axia bindet das wie CalculiX über die Set-Mitgliedschaft.

### Membran (3 translatorische DOF)

| Typ | Knoten | Bemerkung |
|---|---|---|
| M3D4 / M3D4R | 4 | Plane-Stress in der Tangentialebene |
| M3D8 | 8 | quadratisch |
| M3D3 / M3D6 | 3 / 6 | Dreieck |

Dicke über `*MEMBRANE SECTION`. Keine Biegesteifigkeit.

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
| B31 / B21 | 2 | Ende 1, Ende 2 |
| B32 / B22 | 3 | **Ende 1, Ende 2, Mitte** |

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

### Fachwerk und Feder

| Typ | Knoten | Bemerkung |
|---|---|---|
| T3D2 | 2 | 3D-Fachwerk, Querschnitt = erste Zeile `*SOLID SECTION` |
| T3D3 | 3 | quadratischer Stab |
| SPRINGA | 2 | axiale Feder, Steifigkeit über `*SPRING, ELSET=` |
| MASS | 1 | konzentrierte translatorische Masse, `*MASS, ELSET=` |
| ROTARYI | 1 | Drehträgheit I11/I22/I33, `*ROTARY INERTIA, ELSET=` (6 DOF) |
| DASHPOTA | 2 | axialer Dämpfer, `*DASHPOT, ELSET=` (Newmark \(a_1 C\)) |
| GAPUNI | 2 | Spalt: geschlossen (`clearance≤0`) als Axialfeder `*GAP` |

C3D20RI wird als C3D20R gelesen, B21/B22 als B31/B32.

`*EQUATION` koppelt DOFs (Slave-Elimination). `*INCLUDE, INPUT=datei.inp` zieht relative Dateien (CLI, nicht WASM).

`*TIE` (knotenweise, nächster Nachbar; quadratische Flächen inkl. Mittelknoten), `*RIGID BODY, REF NODE=`, `*COUPLING` + `*DISTRIBUTING`/`*KINEMATIC`, `*TRANSFORM, TYPE=R` (lokale DOFs, Ausgabe global).

### Nichtlinear (1.0 T3D2, 1.4 NLGEOM, 1.5 J2)

```
*STEP, NLGEOM
*STATIC
```

Korotationaler Fachwerkstab (T3D2): Materialtangent \(E_t A/L_0\) plus geometrische Steifigkeit \(N/L\,(I-nn^T)\). Kleine Dehnung reproduziert \(u=FL/EA\).

Kontinuum (C3D*): Total-Lagrange, St. Venant–Kirchhoff \(S=\lambda\,\mathrm{tr}(E)\,I+2\mu E\). Kleine Dehnung → linearer Patch (\(u_x=0{,}01\), \(\sigma_{xx}=210\)). Cauchy-Ausgabe \(\sigma=J^{-1}FSF^T\). C3D8I ohne inkompatible Moden (wie C3D8). Gemischte Netze (Schale/Balken+Kontinuum) mit NLGEOM werden abgelehnt.

```
*PLASTIC
210, 0.0
420, 0.01
```

CalculiX-Reihenfolge \(\sigma_y, \bar\varepsilon^p\). T3D2: 1D-J2. Kontinuum: Radial-Return, \(q=\sqrt{3J_2}\), isotrope Verfestigung, konsistente Tangente. Uniaxial C3D8 mit \(\sigma=250\), \(H=21000\) liefert \(u_x\approx 0{,}03095\). FRD-Block `PEEQ`. `NLGEOM`+`*PLASTIC` auf C3D*: J2 auf Green–Lagrange, PK2 plus geometrische Steifigkeit (kleine Dehnung wie ohne NLGEOM).

### Kontakt (1.6)

```
*SURFACE, NAME=MASTER, TYPE=ELEMENT
1, S2
*SURFACE, NAME=SLAVE, TYPE=ELEMENT
2, S1
*SURFACE INTERACTION, NAME=INT
*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR
1e8
*CONTACT PAIR, INTERACTION=INT, TYPE=NODE TO SURFACE
SLAVE, MASTER
```

Penalty, Node-to-Surface. Slave-Knoten gegen Master-Flächen (C3D8/C3D20/C3D4/C3D10/C3D6). Aktiv bei \(g=n\cdot(x_s-x_c)\le 0\). Zwei Würfel in Serie: \(u_{z,\mathrm{oben}}=-0{,}01\) → Interface \(-0{,}005\).

`*FRICTION` — Coulomb auf Knotenkräften, small sliding. Haften wenn \(|F_t|\le\mu |F_n|\), sonst Gleiten. Würfel auf Fundament, \(\mu=0.8\): Scherung bei haftendem Interface. Nicht mit `NLGEOM` oder `*PLASTIC` kombiniert.

### Riks (1.8)

```
*STEP, NLGEOM
*STATIC, RIKS
0.05, 1.0, 1e-4, 0.2, 80
*CLOAD
3, 2, -200
```

Datenzeile wie CalculiX: Anfangsinkrement \(\Delta\lambda\), Period (Ziel-\(\lambda\)), min, max, max. Inkremente. Last \(F=\lambda F_{\mathrm{ref}}\). Constraint: Crisfield zylindrisch \(\Delta u\cdot\Delta u=\Delta\ell^2\), Vorzeichen aus \(v\cdot\Delta u_{\mathrm{prev}}\). Limitpunkt / Snap-Through (von Mises-Fachwerk) wird durchlaufen; \(\lambda\) darf fallen. Ausgabe: \(\lambda\), Inkrementzahl. Nicht mit Kontakt kombiniert.


### Wärmeleitung und Dynamik (1.1)

```
*CONDUCTIVITY
50
*BOUNDARY
1, 11, 11, 0
3, 11, 11, 100
*STEP
*HEAT TRANSFER, STEADY STATE
```

Stationäre Leitung \(K T = F\). NT ist DOF 11. `*DFLUX` (`BF`, `S1`…`S6`), `*FILM` (`F1`…`F6`), `*CFLUX`. Elemente: T3D2, B31, C3D8, CPS4. FRD-Block `NDTEMP`.

```
*INITIAL CONDITIONS, TYPE=DISPLACEMENT
2, 1, 0.01
*AMPLITUDE, NAME=RAMP
0, 0
1, 1
*DAMPING, ALPHA=0, BETA=0
*STEP
*DYNAMIC
0.005, 0.314159265
*CLOAD, AMPLITUDE=RAMP
2, 1, 10
```

Newmark (mittlere Beschleunigung). Lumped mass wie `*FREQUENCY`. Lasten mit `*AMPLITUDE` zeitabhängig. Ausgabe: letzter Zeitschritt.

## Beispiele und Checks

| Beispiel | Erwartung |
|---|---|
| Zugstab C3D8 / C3D8I / C3D20 | Patch-Test: \(u_x = 0{,}01\), \(\sigma_{xx} = 210\) |
| Wedge C3D6 | Patch-Test entlang der Prisma-Achse: \(u_z = 0{,}01\) |
| C3D15 confined | \(\varepsilon_z=0{,}001\), \(\sigma_z \approx 282{,}7\) |
| CAX4-Rohr | Lamé \(u_r(a) \approx 9{,}08\cdot 10^{-4}\) |
| M3D4-Membran | \(u_x=0{,}01\), \(\sigma_{xx}=210\) (\(\nu=0\)) |
| T3D2 / `*EQUATION` | \(u = FL/EA\); zwei Stäbe in Serie \(u_{\mathrm{tip}}=0{,}5\) |
| T3D2 NLGEOM | kleine Dehnung: \(u=0{,}5\) |
| C3D8 NLGEOM Patch | \(u_x=0{,}01\), \(\sigma_{xx}=210\) |
| C3D8 SVK \(\lambda=1{,}2\) | Cauchy \(\sigma_{xx}=26{,}4\) |
| T3D2 `*PLASTIC` | \(\sigma=250\), \(u\approx 3{,}095\) (über elastisch \(1{,}19\)) |
| C3D8 `*PLASTIC` | \(\sigma=250\), \(u_x\approx 0{,}03095\) |
| zwei C3D8 Kontakt | \(u_{z,\mathrm{Interface}}=-0{,}005\) |
| Coulomb-Haften \(\mu=0{,}8\) | \(u_{x,\mathrm{slave}}\approx 0\) |
| C3D8 NLGEOM+J2 | \(u_x\approx 0{,}03095\) |
| T3D2 `*STATIC, RIKS` | Snap-Through, Apex \(u_y < -0{,}8\), \(\lambda\approx 1\) |
| C3D8 Riks Patch | \(u_x=0{,}01\), \(\lambda=1\) |
| Wärmeleitung T3D2 | \(T(0)=0\), \(T(L)=100\) → \(T(L/2)=50\) |
| SDOF `*DYNAMIC` | \(k=100\), \(m=1\), \(u(T/2)=-u_0\) |
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
