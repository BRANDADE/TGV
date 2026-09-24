"""A3 — ce que TGV ne sait pas classer n'est ni 'normal' ni supprimé de l'export."""
from trgt_vcf import Call, record

from scripts.bio.labels import NOTE_MOTIF_DISCORDANCE, NOTE_OUT_OF_RANGE, UNCLASSIFIED
from scripts.ui.html_export import generate_html_table

RFC1_MOTIFS = ["AAAAG", "AAAGG", "AAGGG", "ACAGG", "AGGGC"]


def test_value_outside_all_ranges_is_unclassified(analyze):
    # SCA36 : 15–649 absent du YAML (GeneReviews : signification incertaine)
    line = record("SCA36_NOP56", ["GGCCTG"], [Call("GGCCTG" * 8), Call("GGCCTG" * 100)])
    r = analyze([line])["SCA36_NOP56"]
    assert r.classification1_raw == "normal"
    assert r.classification2_raw == UNCLASSIFIED
    assert r.classification2_note == NOTE_OUT_OF_RANGE


def test_overlapping_ranges_are_unclassified(analyze):
    # FMR1 : 200 est à la fois 'premutation' [55, 200] et 'pathogenic_complete' [200, null]
    line = record("FXS_FMR1", ["CGG"], [Call("CGG" * 30), Call("CGG" * 200)])
    r = analyze([line])["FXS_FMR1"]
    assert r.classification2_raw == UNCLASSIFIED
    assert "premutation" in r.classification2_note and "pathogenic_complete" in r.classification2_note


def test_ranges_around_overlap_still_classified(analyze):
    line = record("FXS_FMR1", ["CGG"], [Call("CGG" * 199), Call("CGG" * 201)])
    r = analyze([line])["FXS_FMR1"]
    assert (r.classification1_raw, r.classification2_raw) == ("premutation", "pathogenic_complete")


def test_motif_discordance_keeps_locus_as_unclassified(analyze):
    # Motifs du catalogue sans rapport avec les groupes du YAML (CAG/CAA)
    line = record("SCA17_TBP", ["AAT"], [Call("AAT" * 20), Call("AAT" * 30)])
    r = analyze([line])["SCA17_TBP"]
    assert r.has_clinical
    assert (r.classification1_raw, r.classification2_raw) == (UNCLASSIFIED, UNCLASSIFIED)
    assert r.classification1_note == NOTE_MOTIF_DISCORDANCE
    assert r.display_row.genotype == ". / ."


def test_rfc1_without_normal_range_shows_other_motif(analyze):
    # RFC1 : pas de plage 'normal' dans le YAML → AAAAG(11) 'unclassified' ;
    # le génotype affiche toujours le motif majoritaire (full_with_others).
    line = record("CANVAS_RFC1", RFC1_MOTIFS, [Call("AAAAG" * 11), Call("AAGGG" * 450)])
    r = analyze([line])["CANVAS_RFC1"]
    assert (r.classification1_raw, r.classification2_raw) == (UNCLASSIFIED, "pathogenic")
    assert (r.genotype1_raw, r.genotype2_raw) == ("11 (AAAAG)", "450 (AAGGG)")


def test_unclassified_loci_are_exported_with_explanation(analyze):
    lines = [
        record("SCA36_NOP56", ["GGCCTG"], [Call("GGCCTG" * 8), Call("GGCCTG" * 100)], pos=1000),
        record("SCA17_TBP", ["AAT"], [Call("AAT" * 20), Call("AAT" * 30)], pos=5000),
    ]
    results = analyze(lines)
    rows = [
        {"Locus": r.locus, "Profondeur": r.display_row.depth, "Génotype": r.display_row.genotype,
         "Classification": r.display_row.classification, "Result_obj": r,
         "Classification_auto": r.display_row.classification, "Genotype_auto": r.display_row.genotype}
        for r in results.values()
    ]
    html = generate_html_table(["Locus", "Profondeur", "Génotype", "Classification"], rows, "S1", run_id="R1")
    assert "SCA36 (NOP56)" in html and "SCA17 (TBP)" in html
    assert "Allèle 2 non classé : valeur hors des plages du YAML" in html
    assert "Allèle 1 non classé : motifs TRGT absents" in html
