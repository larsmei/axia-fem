# NAFEMS-Benchmarks

Fälle aus *The Standard NAFEMS Benchmarks* (TNSB, Rev. 3) und der Wärmeplatte
von Cameron, Casey & Simpson. Decks erzeugt `solver/scripts/gen_nafems.py`.
Tests: `cargo test --features cli --lib nafems`.

| Test | Elemente | Datei | Ziel | Axia |
|---|---|---|---|---|
| LE1 elliptische Membran | CPS8 | `le1_cps8.inp` | σyy(D) = 92.7 MPa | 93.39 MPa (+0.7 %), ux(D) = −0.102 mm |
| LE1, Dreiecke | CPS6 | `le1_cps6.inp` | σyy(D) = 92.7 MPa | 89.40 MPa (−3.6 %) |
| LE6 schiefe Platte | S8 | `le6_s8.inp` | σ1(E) = 0.802 MPa | 0.771 MPa (−3.9 %) |
| LE10 dicke Platte | C3D20 | `le10_c3d20.inp` | σyy(D) = −5.38 MPa | −5.51 MPa (+2.4 %) |
| LE10, verfeinert | C3D20 | `le10_c3d20_fine.inp` | σyy(D) = −5.38 MPa | −5.50 MPa (+2.3 %) |
| T4 Konvektion | C3D8 | `t4_c3d8.inp` | T(0.6 m, 0.2 m) = 18.3 °C | 18.21 °C |

LE1: Viertelsymmetrie, Außenkante 10 MPa nach außen, D = innerer Punkt auf +x. Längen in mm.
LE6: Parallelogramm, Seiten 1 m, Schiefe 30°, t = 10 mm, p = −0.7 kPa, gelenkig gelagert. Die Schale speichert die Faser z = +t/2; die Zug-Hauptspannung der Unterseite ist der Betrag der Druckfaser. E = Plattenmitte.
LE10: Grundfläche wie LE1, Dicke 600 mm, 1 MPa auf der Oberseite. D = oberer innerer Punkt.
T4: 0.6×1.0 m, k = 52, Boden 100 °C, oben/rechts Film h = 750 gegen 0 °C. Dünne C3D8-Scheibe, weil `*FILM` an Hex-Flächen hängt.

Nicht in der Suite, obwohl spezifiziert: LE3 (Halbkugelschale, Ziel ux(A) = 0.185 m) — S4R bleibt bei 24×24 bei 0.070 m (Membran-Locking). LE5 (Z-Querschnitt, Schale/Balken, −108 MPa) und LE11 (Temperaturspannung, −105 MPa) brauchen Querschnittsdaten bzw. Wärmeausdehnung auf C3D20, die der Solver so nicht hat. LE7/LE8 sind Rotationsschalen, kein CAX-Kontinuum.
