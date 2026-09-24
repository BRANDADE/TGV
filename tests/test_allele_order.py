"""A4 — allèles ordonnés en amont (plus petit génotype clinique en premier) ; plus de bouton ⇅."""
from trgt_vcf import Call, rc, record

from scripts.ui.html_export import generate_html_table

RFC1_MOTIFS = ["AAAAG", "AAAGG", "AAGGG", "ACAGG", "AGGGC"]


def test_larger_allele_first_in_vcf_is_displayed_second(analyze):
    # Ordre TRGT : 42 (référence) puis 30 ; l'affichage doit donner 30 / 42.
    line = record(
        "SCA1_ATXN1", ["CTG", "TTG"],
        [Call(rc("CAG" * 42), depth=27), Call(rc("CAG" * 14 + "CAT" + "CAG" * 15), depth=31)],
    )
    r = analyze([line])["SCA1_ATXN1"]
    assert r.display_row.genotype == "30 / 42"
    assert r.display_row.classification == "normal / pathogenic_complete"
    # Les champs suivent l'allèle : 31 lectures pour le 30, 27 pour le 42
    assert r.display_row.depth == "⚠ 31 / ⚠ 27"
    assert r.display_row.size == "90 / 126"


def test_rfc1_order_uses_displayed_genotype(analyze):
    line = record("CANVAS_RFC1", RFC1_MOTIFS, [Call("AAGGG" * 450), Call("AAAAG" * 11)])
    r = analyze([line])["CANVAS_RFC1"]
    assert r.display_export.genotype == "11 (AAAAG) / 450 (AAGGG)"
    assert r.display_row.classification == "unclassified / pathogenic"


def test_pure_only_orders_on_first_displayed_number(analyze):
    # FMR1 (pure_only) : '31 (32)' vs '30 (30)' → 30 d'abord
    line = record(
        "FXS_FMR1", ["CGG", "AGG"],
        [Call("CGG" * 10 + "AGG" + "CGG" * 21), Call("CGG" * 30)],
    )
    r = analyze([line])["FXS_FMR1"]
    assert r.display_row.genotype == "30 (30) / 31 (32)"


def test_haploid_called_allele_always_first(analyze):
    line = record("FXS_FMR1", ["CGG"], [Call("CGG" * 60)], ref="CGG" * 20, gt="1")
    r = analyze([line])["FXS_FMR1"]
    assert r.display_row.genotype == "60 (60) / -"


def test_non_clinical_locus_ordered_by_length(analyze):
    line = record("HD_HTT", ["CAG"], [Call("CAG" * 40), Call("CAG" * 17)])
    r = analyze([line])["HD_HTT"]
    assert r.display_row.size == "51 / 120"


def test_html_export_has_no_swap_button(analyze):
    line = record("SCA17_TBP", ["CAG", "CAA"], [Call("CAG" * 36), Call("CAG" * 38)])
    r = analyze([line])["SCA17_TBP"]
    row = {"Locus": r.locus, "Profondeur": r.display_row.depth, "Génotype": r.display_row.genotype,
           "Classification": r.display_row.classification, "Result_obj": r,
           "Classification_auto": r.display_row.classification, "Genotype_auto": r.display_row.genotype}
    html = generate_html_table(["Locus", "Profondeur", "Génotype", "Classification"], [row], "S1")
    assert "swapRowAlleles" not in html and "swapCompositeValue" not in html and "⇅" not in html
    assert ".slice(1).map" in html
