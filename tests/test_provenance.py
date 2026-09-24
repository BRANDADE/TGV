"""E1/E2 — traçabilité : versions, catalogue, empreinte du YAML, surcharges manuelles."""
import csv
import re

from trgt_vcf import Call, record, vcf_text, write_zip

import tgv_cli
from scripts.core.comments import result_comments
from scripts.core.provenance import catalog_from_command, format_source, read_vcf_header, tgv_provenance
from scripts.core.version import build_commit
from scripts.ui.html_export import generate_html_table


def _run(tmp_path, version="5.1.0"):
    line = record("SCA17_TBP", ["CAG", "CAA"], [Call("CAG" * 36), Call("CAG" * 52)])
    return write_zip(tmp_path / "RUN9-trgt_vcfs.zip", {
        "S1.trgt.sorted.vcf": vcf_text([line], version=version, repeats="/data/bed/panel_v3.bed"),
    })


def test_read_vcf_header(tmp_path):
    info = read_vcf_header(_run(tmp_path), "S1.trgt.sorted.vcf")
    assert info["trgt_version"] == "5.1.0"
    assert info["catalog"] == "panel_v3.bed"
    assert "--preset targeted" in info["trgt_command"]


def test_catalog_from_command():
    assert catalog_from_command("trgt genotype --genome g.fa --repeats /x/y.bed --reads r.bam") == "y.bed"
    assert catalog_from_command("trgt genotype -b cat.bed") == "cat.bed"
    assert catalog_from_command("") == ""


def test_format_source(thresholds_data):
    assert format_source(thresholds_data["SCA17_TBP"]["source"]) == "HAS - Volet 1 (02-2025)"
    assert format_source(thresholds_data["FRDA_FXN"]["source"]) == "HAS - Volet 1 (02-2025); Benkirane et al. (2025)"


def test_tgv_provenance():
    p = tgv_provenance()
    assert p["tgv_version"]
    assert build_commit()
    assert re.fullmatch(r"[0-9a-f]{64}", p["thresholds_sha256"])


def test_cli_rows_carry_provenance(tmp_path):
    out = str(tmp_path / "export.tsv")
    assert tgv_cli.main(["--zip", _run(tmp_path), "--out", out]) == 0
    with open(out, newline="", encoding="utf-8") as f:
        row = next(csv.DictReader(f, delimiter="\t"))
    assert row["trgt_version"] == "5.1.0"
    assert row["catalog"] == "panel_v3.bed"
    assert row["thresholds_source"] == "HAS - Volet 1 (02-2025)"
    assert row["tgv_version"] and row["tgv_commit"]
    assert re.fullmatch(r"[0-9a-f]{64}", row["thresholds_sha256"])


def test_cli_rejects_mismatching_trgt_version(tmp_path):
    out = str(tmp_path / "export.tsv")
    assert tgv_cli.main(["--zip", _run(tmp_path), "--out", out, "--trgt-version", "4.0.0"]) == tgv_cli.EXIT_FAILURE


def test_manual_overrides_are_traced_in_comments(analyze):
    line = record("SCA17_TBP", ["CAG", "CAA"], [Call("CAG" * 36), Call("CAG" * 52)])
    r = analyze([line])["SCA17_TBP"]
    assert result_comments(r) == []

    r.classification2_bio = "pathogenic_incomplete"
    r.classification1_bio = "normal"
    r.genotype2_bio = 51
    comments = result_comments(r)
    assert "Classification modifiée manuellement (auto : normal / pathogenic_complete)" in comments
    assert "Génotype modifié manuellement (auto : 36 / 52)" in comments

    # Valeurs identiques à l'automatique : pas de trace
    r.classification2_bio = "pathogenic_complete"
    r.genotype1_bio, r.genotype2_bio = 36, 52
    assert result_comments(r) == []


def test_html_export_footer_shows_provenance(analyze):
    line = record("SCA17_TBP", ["CAG", "CAA"], [Call("CAG" * 36), Call("CAG" * 52)])
    r = analyze([line])["SCA17_TBP"]
    row = {"Locus": r.locus, "Profondeur": r.display_row.depth, "Génotype": r.display_row.genotype,
           "Classification": r.display_row.classification, "Result_obj": r,
           "Classification_auto": r.display_row.classification, "Genotype_auto": r.display_row.genotype}
    provenance = {**tgv_provenance(), "trgt_version": "5.1.0", "catalog": "panel_v3.bed"}
    page = generate_html_table(["Locus"], [row], "S1", provenance=provenance)
    assert "TRGT 5.1.0" in page and "catalogue panel_v3.bed" in page
    assert provenance["thresholds_sha256"][:16] in page
