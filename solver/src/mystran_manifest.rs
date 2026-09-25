//! MYSTRAN bulk-card contract.
//!
//! The list is the User Manual (2025-09-22, chapter 11) crossed with the
//! decks under `Build_Test_Cases`. A card is implemented, declined with a
//! hard error, or a known MYSTRAN defect that Axia does not reproduce.
//! Nothing in the manual stays planned, and a bulk name outside the manual
//! is a hard error.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CardStatus {
    Implemented,
    Planned,
    Declined,
    /// Present so the status exists. Axia solves the case correctly instead
    /// of copying the MYSTRAN defect. See the card note.
    MystranBug,
}

pub struct Card {
    pub name: &'static str,
    pub status: CardStatus,
    pub note: &'static str,
}

pub static CARDS: &[Card] = &[
    card("PARAM", CardStatus::Implemented, "WTMASS, AUTOSPC, K6ROT, GRDPNT"),
    card("DEBUG", CardStatus::Implemented, "akzeptiert, ohne Wirkung"),
    card("EIGRL", CardStatus::Implemented, "alle Eigenvektoren in der F06"),
    card("GRDSET", CardStatus::Implemented, ""),
    card("GRID", CardStatus::Implemented, ""),
    card("CORD2R", CardStatus::Implemented, ""),
    card("CORD2C", CardStatus::Implemented, "R, θ in Grad, Z"),
    card("CORD2S", CardStatus::Implemented, "R, θ, φ in Grad"),
    card("CORD1R", CardStatus::Implemented, ""),
    card("CORD1C", CardStatus::Implemented, ""),
    card("CORD1S", CardStatus::Implemented, ""),
    card("MAT1", CardStatus::Implemented, ""),
    card("PSHELL", CardStatus::Implemented, "12I/T³ skaliert die Biegesteifigkeit"),
    card("PSOLID", CardStatus::Implemented, ""),
    card("PROD", CardStatus::Implemented, ""),
    card("PBAR", CardStatus::Implemented, ""),
    card("PBARL", CardStatus::Implemented, "ROD, TUBE, TUBE2, BAR, BOX, I, T, L"),
    card("CROD", CardStatus::Implemented, "F06 Normalkraft und Axialspannung, Zug positiv"),
    card("CONROD", CardStatus::Implemented, "wie CROD"),
    card("CBAR", CardStatus::Implemented, "Hermite + Schub, PA/PB, Versatz, F06-Schnittgrößen"),
    card("CBEAM", CardStatus::Implemented, "wie CBAR"),
    card("CQUAD4", CardStatus::Implemented, "MITC4, THETA und MCID"),
    card("CQUAD4K", CardStatus::Implemented, "wie CQUAD4"),
    card("CTRIA3", CardStatus::Implemented, "THETA und MCID"),
    card("CTRIA3K", CardStatus::Implemented, "wie CTRIA3"),
    card("CTETRA", CardStatus::Implemented, ""),
    card("CHEXA", CardStatus::Implemented, ""),
    card("CPENTA", CardStatus::Implemented, ""),
    card("CELAS1", CardStatus::Implemented, "über PELAS, komponentenweise"),
    card("CELAS2", CardStatus::Implemented, "Komponente, auch gegen Erde"),
    card("CELAS3", CardStatus::Implemented, "skalar über PELAS"),
    card("CELAS4", CardStatus::Implemented, "skalar"),
    card("PELAS", CardStatus::Implemented, ""),
    card("CMASS1", CardStatus::Implemented, "über PMASS"),
    card("CMASS2", CardStatus::Implemented, ""),
    card("CMASS3", CardStatus::Implemented, "skalar über PMASS"),
    card("CMASS4", CardStatus::Implemented, "skalar"),
    card("PMASS", CardStatus::Implemented, ""),
    card("CONM2", CardStatus::Implemented, "Versatz und Drehträgheit"),
    card("CSHEAR", CardStatus::Implemented, "nur Schub"),
    card("PSHEAR", CardStatus::Implemented, ""),
    card("CBUSH", CardStatus::Implemented, ""),
    card("PBUSH", CardStatus::Implemented, "Steifigkeit K"),
    card("BAROR", CardStatus::Implemented, "Vorgabe für folgende CBAR"),
    card("RBE2", CardStatus::Implemented, "CM teilweise"),
    card("RBE3", CardStatus::Implemented, "gewichtete Interpolation"),
    card("MPC", CardStatus::Implemented, ""),
    card("MPCADD", CardStatus::Implemented, ""),
    card("FORCE", CardStatus::Implemented, ""),
    card("MOMENT", CardStatus::Implemented, ""),
    card("PLOAD2", CardStatus::Implemented, "positiv entgegen der Normalen"),
    card(
        "PLOAD4",
        CardStatus::MystranBug,
        "Axia wertet PLOAD4 auch in LOAD aus. MYSTRAN Issue 196 liefert dort Null.",
    ),
    card("GRAV", CardStatus::Implemented, ""),
    card("LOAD", CardStatus::Implemented, ""),
    card("RFORCE", CardStatus::Implemented, "ω = 2π·V Umdrehungen"),
    card("TEMP", CardStatus::Implemented, "Gitter, überschreibt TEMPD"),
    card("TEMPD", CardStatus::Implemented, "Vorgabe für alle Gitter"),
    card("TEMPP1", CardStatus::Implemented, "mittlere Plattensemperatur; TPRIME abgewiesen"),
    card("TEMPRB", CardStatus::Implemented, "mittlere Stabtemperatur; Gradienten abgewiesen"),
    card("SPC", CardStatus::Implemented, "auch erzwungene Verschiebung"),
    card("SPC1", CardStatus::Implemented, ""),
    card("SPCADD", CardStatus::Implemented, ""),
    card("ASET", CardStatus::Declined, "ASET-Reduktion wird nicht still verworfen. Das Deck wird abgewiesen."),
    card("ASET1", CardStatus::Declined, "ASET1-Reduktion wird abgewiesen."),
    card("OMIT", CardStatus::Declined, "OMIT-Reduktion wird abgewiesen."),
    card("OMIT1", CardStatus::Declined, "OMIT1-Reduktion wird abgewiesen."),
    card("EIGR", CardStatus::Implemented, "Subspace statt Givens; NORM MAX, MASS, POINT"),
    card("MAT2", CardStatus::Implemented, "anisotrope Scheibe"),
    card("MAT8", CardStatus::Implemented, "orthotrop, Winkel wird angesetzt"),
    card("MAT9", CardStatus::Implemented, "nur lineares CHEXA"),
    card(
        "PCOMP",
        CardStatus::Implemented,
        "lineare Steifigkeit, Versagenskriterien ohne Wirkung",
    ),
    card("PCOMP1", CardStatus::Implemented, "wie PCOMP, gleiche Lagen"),
    card(
        "CUSERIN",
        CardStatus::Declined,
        "Superelement CUSERIN liegt außerhalb des linearen Elementkerns.",
    ),
    card(
        "PARVEC",
        CardStatus::Declined,
        "PARVEC ist ein MYSTRAN-Debug-Vektor und wird nicht gerechnet.",
    ),
    card(
        "PARVEC1",
        CardStatus::Declined,
        "PARVEC1 ist ein MYSTRAN-Debug-Vektor und wird nicht gerechnet.",
    ),
    card("CQUAD8", CardStatus::Implemented, "S8, Ecken dann Mitten, THETA/MCID"),
    card("PLOTEL", CardStatus::Implemented, "nur Plot, ohne Steifigkeit"),
    card(
        "PUSERIN",
        CardStatus::Declined,
        "Superelement PUSERIN liegt außerhalb des linearen Elementkerns.",
    ),
    card("RSPLINE", CardStatus::Declined, "RSPLINE wird nicht gerechnet."),
    card("SEQGP", CardStatus::Implemented, "Sequenzierung ohne Wirkung"),
    card("SLOAD", CardStatus::Declined, "skalare Last wird nicht gerechnet."),
    card("SPOINT", CardStatus::Implemented, "skalarer Punkt für CELAS3/CMASS3"),
    card(
        "SUPORT",
        CardStatus::Declined,
        "SUPORT (Craig-Bampton) wird nicht gerechnet.",
    ),
    card("USET", CardStatus::Declined, "USET-Menge wird abgewiesen."),
    card("USET1", CardStatus::Declined, "USET1-Menge wird abgewiesen."),
];

const fn card(name: &'static str, status: CardStatus, note: &'static str) -> Card {
    Card { name, status, note }
}

pub fn lookup(name: &str) -> Option<&'static Card> {
    CARDS.iter().find(|c| c.name.eq_ignore_ascii_case(name))
}

pub fn names_with(status: CardStatus) -> Vec<&'static str> {
    let mut v: Vec<_> = CARDS
        .iter()
        .filter(|c| c.status == status)
        .map(|c| c.name)
        .collect();
    v.sort_unstable();
    v
}

/// Cards that are neither implemented, declined, nor a documented MYSTRAN defect.
pub fn open_gaps() -> Vec<&'static str> {
    names_with(CardStatus::Planned)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn manifest_is_unique_and_classified() {
        let mut names: Vec<_> = CARDS.iter().map(|c| c.name).collect();
        let n = names.len();
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), n, "doppelte Kartennamen");
        assert!(CARDS.iter().any(|c| c.status == CardStatus::MystranBug));
        assert!(CARDS.iter().any(|c| c.status == CardStatus::Declined));
        assert!(!names_with(CardStatus::Implemented).is_empty());
    }

    #[test]
    fn mystran_manifest_closed() {
        let gaps = open_gaps();
        assert!(gaps.is_empty(), "noch offen: {}", gaps.join(", "));
    }

    /// Chapter 11 of the 2025-09-22 manual. Every name is classified.
    #[test]
    fn manual_2025_09_22_bulk_is_classified() {
        const MANUAL: &[&str] = &[
            "ASET", "ASET1", "BAROR", "CBAR", "CBUSH", "CELAS1", "CELAS2", "CELAS3", "CELAS4",
            "CHEXA", "CMASS1", "CMASS2", "CMASS3", "CMASS4", "CONM2", "CONROD", "CORD1C",
            "CORD1R", "CORD1S", "CORD2C", "CORD2R", "CORD2S", "CPENTA", "CQUAD4", "CQUAD4K",
            "CQUAD8", "CROD", "CSHEAR", "CTETRA", "CTRIA3", "CTRIA3K", "CUSERIN", "DEBUG", "EIGR",
            "EIGRL", "FORCE", "GRAV", "GRDSET", "GRID", "LOAD", "MAT1", "MAT2", "MAT8", "MAT9",
            "MOMENT", "MPC", "MPCADD", "OMIT", "OMIT1", "PARAM", "PARVEC", "PARVEC1", "PBAR",
            "PBARL", "PBUSH", "PCOMP", "PCOMP1", "PELAS", "PLOAD2", "PLOAD4", "PLOTEL", "PROD",
            "PSHEAR", "PSHELL", "PSOLID", "PUSERIN", "RBE2", "RBE3", "RFORCE", "RSPLINE", "SEQGP",
            "SLOAD", "SPC", "SPC1", "SPCADD", "SPOINT", "SUPORT", "TEMP", "TEMPD", "TEMPP1",
            "TEMPRB", "USET", "USET1",
        ];
        for name in MANUAL {
            let card = lookup(name).unwrap_or_else(|| panic!("{name} fehlt im Manifest"));
            assert_ne!(card.status, CardStatus::Planned, "{name} ist noch geplant");
        }
    }
}
