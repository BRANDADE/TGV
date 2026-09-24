"""A2 (GT phasé) et B1 (no-call, haploïde, MS vide) sur des sorties au format TRGT 5.x."""
import pytest
from trgt_vcf import Call, rc, record

from scripts.core.vcf_parser import parse_gt

# SCA1 (orientation clinique RC) : catalogue du CHU au brin + (CTG/TTG).
SCA1_MOTIFS = ["CTG", "TTG"]
SCA1_INTERRUPTED_30 = "CAG" * 14 + "CAT" + "CAG" * 15   # 29 CAG + 1 CAT → 30
SCA1_PURE_42 = "CAG" * 42


@pytest.mark.parametrize("gt, expected", [
    ("0/1", ([0, 1], 2)),
    ("1|0", ([0, 1], 2)),
    ("2|1", ([1, 2], 2)),
    ("1/1", ([1, 1], 2)),
    ("1", ([1], 1)),
    (".", ([], 0)),
    ("", ([], 0)),
])
def test_parse_gt_returns_trgt_allele_order(gt, expected):
    assert parse_gt(gt) == expected


@pytest.mark.parametrize("gt", ["0/1", "0|1", "1|0"])
def test_sca1_phased_genotype_keeps_sequence_format_pairing(analyze, gt):
    # Réf. revue A2 : écrit '1|0', l'ancien code donnait 29/43 normal/normal (faux négatif).
    line = record(
        "SCA1_ATXN1", SCA1_MOTIFS,
        [Call(rc(SCA1_INTERRUPTED_30), depth=31), Call(rc(SCA1_PURE_42), depth=27)],
        gt=gt,
    )
    r = analyze([line])["SCA1_ATXN1"]
    assert (r.genotype1_raw, r.genotype2_raw) == ("30", "42")
    assert (r.classification1_raw, r.classification2_raw) == ("normal", "pathogenic_complete")
    assert (r.depth1_raw, r.depth2_raw) == ("31", "27")


def test_no_call_does_not_crash(analyze):
    line = record("SCA2_ATXN2", ["CTG", "TTG"], [], ref=rc("CAG" * 22))
    r = analyze([line])["SCA2_ATXN2"]
    assert r.display_row.classification == "no_call / no_call"
    assert r.display_row.genotype == "no_call / no_call"
    assert r.display_row.depth == ". / ."
    assert r.display_export.depth1 == "."


def test_haploid_allele_does_not_crash(analyze):
    # Format TRGT réel (caryotype XY) : GT '1', un seul jeu de valeurs FORMAT.
    line = record("FXS_FMR1", ["CGG"], [Call("CGG" * 11, depth=29)], ref="CGG" * 30, gt="1")
    assert line.split("\t")[9].startswith("1:33:33-33:29:11:0(0-33):1:")

    r = analyze([line])["FXS_FMR1"]
    assert r.status2 == "absent"
    assert r.display_row.classification == "normal / -"
    assert r.display_row.genotype == "11 (11) / -"
    assert r.display_row.depth == "⚠ 29 / -"
    assert r.display_row.size == "33 / -"


def test_empty_motif_segmentation_does_not_crash(analyze):
    # MS='.' pour un allèle appelé (FW) : la classification repose sur MC.
    line = record(
        "SCA17_TBP", ["CAG", "CAA"],
        [Call("CAG" * 38, ms="."), Call("CAG" * 52)],
    )
    r = analyze([line])["SCA17_TBP"]
    assert (r.genotype1_raw, r.genotype2_raw) == ("38", "52")
    assert (r.classification1_raw, r.classification2_raw) == ("normal", "pathogenic_complete")
    assert r.seg1_raw == ""


def test_low_depth_marker_only_in_display(analyze):
    line = record("SCA17_TBP", ["CAG", "CAA"], [Call("CAG" * 36, depth=20), Call("CAG" * 38, depth=80)])
    r = analyze([line])["SCA17_TBP"]
    assert r.display_row.depth == "⚠ 20 / 80"
    assert (r.display_export.depth1, r.display_export.depth2) == ("20", "80")
