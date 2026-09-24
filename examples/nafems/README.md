# NAFEMS-Benchmarks

Drei Fälle aus *The Standard NAFEMS Benchmarks* (TNSB, Rev. 3) bzw. der
Wärmeplatte von Cameron, Casey & Simpson. Decks erzeugt
`solver/scripts/gen_nafems.py`. Tests: `cargo test --features cli --lib nafems`.

Einheiten LE1/LE10: mm, N, MPa. T4: m, °C, W/(m·K).

| Test | Datei | Netz | Ziel | Axia |
|---|---|---|---|---|
| LE1 elliptische Membran | `le1_cps8.inp` | CPS8, 8×24 | σyy(D) = 92.7 MPa | 93.39 MPa (+0.7 %), ux(D) = −0.102 mm |
| LE10 dicke Platte | `le10_c3d20.inp` | C3D20, 6×4×2 | σyy(D) = −5.38 MPa | −5.51 MPa (+2.4 %) |
| LE10, verfeinert | `le10_c3d20_fine.inp` | C3D20, 12×8×4 | σyy(D) = −5.38 MPa | −5.50 MPa (+2.3 %) |
| T4 Konvektion | `t4_c3d8.inp` | C3D8, 24×40×1 | T(0.6 m, 0.2 m) = 18.3 °C | 18.21 °C (−0.09 °C) |

LE1: Viertelsymmetrie, Außenkante 10 MPa nach außen, D = innerer Punkt auf +x.
LE10: dieselbe Grundfläche, Dicke 600 mm, 1 MPa Druck auf der Oberseite, Außenkante in der Ebene gehalten, z nur in der Mittelebene der Außenkante. D = oberer innerer Punkt.
T4: Platte 0.6×1.0 m, k = 52, Boden 100 °C, links isoliert, oben und rechts Film h = 750 gegen 0 °C. Als dünne C3D8-Scheibe, weil `*FILM` an Hex-Flächen hängt.
