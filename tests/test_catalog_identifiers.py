"""B2 — catalogues publics : TRID sans « _ », alias vers le YAML, cadre de lecture des motifs."""
import os

import yaml
from trgt_vcf import Call, record

from scripts.core.analysis import resolve_panel
from scripts.core.sequence_utils import motif_frame_map

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_aliases():
    with open(os.path.join(ROOT, "configs", "trid_aliases.yaml"), encoding="utf-8") as f:
        return yaml.safe_load(f)


def test_motif_frame_map():
    # ATXN1 (catalogue public) : TGC → RC 'GCA', rotation de 'CAG'
    assert motif_frame_map(["GCA"], {"CAG", "CAA"}) == {"GCA": "CAG"}
    # Motif déjà présent tel quel : rien à renommer
    assert motif_frame_map(["CAG", "CAA"], {"CAG", "CAA"}) == {}
    # Motif sans rapport : pas de renommage
    assert motif_frame_map(["GCC"], {"CAG"}) == {}
    # Le catalogue déclare déjà CAG : GCA reste un motif distinct
    assert motif_frame_map(["CAG", "GCA"], {"CAG"}) == {}


def test_aliases_target_existing_yaml_blocks(thresholds_data):
    for catalog_id, block in load_aliases().items():
        assert isinstance(thresholds_data.get(block), dict), f"{catalog_id} → {block}"


def test_trid_without_underscore_does_not_crash(analyze):
    # Avant : trid_id.split('_', 1) → ValueError sur 'HTT' (mode --all-loci, « Tout cocher »)
    line = record("HTT", ["CAG", "CCG"], [Call("CAG" * 17), Call("CAG" * 40)])
    r = analyze([line])["HTT"]
    assert r.locus == "HTT"
    assert not r.has_clinical


def _public(analyze, trid, motifs, calls, **kwargs):
    thresholds = yaml.safe_load(open(os.path.join(ROOT, "configs", "clinical_thresholds.yaml"), encoding="utf-8"))
    return analyze([record(trid, motifs, calls, **kwargs)], thresholds=thresholds, aliases=load_aliases())[trid]


def test_public_catalog_atxn1_rc_rotation(analyze):
    # Catalogue public : ID=ATXN1;MOTIFS=TGC (brin +). Clinique : CAG, orientation RC.
    r = _public(analyze, "ATXN1", ["TGC"], [Call("TGC" * 30), Call("TGC" * 42)])
    assert r.locus == "SCA1 (ATXN1)"
    assert r.has_clinical
    assert (r.genotype1_raw, r.genotype2_raw) == ("30", "42")
    assert (r.classification1_raw, r.classification2_raw) == ("normal", "pathogenic_complete")


def test_public_catalog_tbp_fw_rotation(analyze):
    # Catalogue public : ID=TBP;MOTIFS=GCA (FW). Clinique : CAG/CAA, orientation FW.
    r = _public(analyze, "TBP", ["GCA"], [Call("GCA" * 38), Call("GCA" * 52)])
    assert (r.genotype1_raw, r.genotype2_raw) == ("38", "52")
    assert (r.classification1_raw, r.classification2_raw) == ("normal", "pathogenic_complete")


def test_public_catalog_atxn2_rc_rotation(analyze):
    # Catalogue public : ID=ATXN2;MOTIFS=GCT → RC 'AGC', rotation de 'CAG'
    r = _public(analyze, "ATXN2", ["GCT"], [Call("GCT" * 22), Call("GCT" * 40)])
    assert (r.genotype1_raw, r.genotype2_raw) == ("22", "40")
    assert (r.classification1_raw, r.classification2_raw) == ("normal", "pathogenic_complete")


def test_public_catalog_exact_motif_after_rc(analyze):
    # C9ORF72 : GGCCCC → RC 'GGGGCC', identique au YAML (alias seul, pas de renommage)
    r = _public(analyze, "C9ORF72", ["GGCCCC"], [Call("GGCCCC" * 2), Call("GGCCCC" * 40)])
    assert r.locus == "FTDALS1 (C9orf72)"
    assert (r.classification1_raw, r.classification2_raw) == ("normal", "pathogenic_complete")


def test_resolve_panel_with_aliases():
    class T:
        def __init__(self, key):
            self.clinical_key = key
    run_trids = {"ATXN1": T("SCA1_ATXN1"), "HTT": T(None), "SCA2_ATXN2": T("SCA2_ATXN2")}
    found, missing = resolve_panel(["SCA1_ATXN1", "SCA2_ATXN2", "SCA3_ATXN3"], run_trids)
    assert found == ["ATXN1", "SCA2_ATXN2"]
    assert missing == ["SCA3_ATXN3"]
