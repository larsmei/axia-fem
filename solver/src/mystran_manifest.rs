//! MYSTRAN bulk-card contract.
//!
//! The list is the User Manual (2025-09-22, chapter 7) crossed with the 2011
//! manual contents and the decks under `Build_Test_Cases`. A card is either
//! implemented, still planned, declined with a hard error, or a known MYSTRAN
//! defect that Axia does not reproduce.

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
    card("PARAM", CardStatus::Implemented, "WTMASS"),
    card("DEBUG", CardStatus::Implemented, "akzeptiert, ohne Wirkung"),
    card("EIGRL", CardStatus::Implemented, "ND Eigenwerte"),
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
    card("CROD", CardStatus::Implemented, ""),
    card("CONROD", CardStatus::Implemented, ""),
    card("CBAR", CardStatus::Implemented, "Hermite + Schub, PA/PB, Versatz"),
    card("CBEAM", CardStatus::Implemented, "wie CBAR"),
    card("CQUAD4", CardStatus::Implemented, "MITC4, Winkel noch ignoriert"),
    card("CQUAD4K", CardStatus::Implemented, "wie CQUAD4"),
    card("CTRIA3", CardStatus::Implemented, ""),
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
    card("EIGR", CardStatus::Planned, ""),
    card("MAT2", CardStatus::Planned, ""),
    card("MAT8", CardStatus::Planned, ""),
    card("MAT9", CardStatus::Planned, ""),
    card("PCOMP", CardStatus::Planned, ""),
    card("PCOMP1", CardStatus::Planned, ""),
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
    #[ignore = "Phase 8: keine geplanten Karten mehr"]
    fn mystran_manifest_closed() {
        let gaps = open_gaps();
        assert!(gaps.is_empty(), "noch offen: {}", gaps.join(", "));
    }
}
