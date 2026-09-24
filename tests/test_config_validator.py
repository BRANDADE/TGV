"""E3 — validation de clinical_thresholds.yaml."""
import copy

from scripts.bio.clinical_config_validator import coverage_gaps, validate_thresholds


def test_shipped_yaml_has_no_error_and_expected_warnings(thresholds_data):
    errors, warnings = validate_thresholds(thresholds_data)
    assert errors == []
    expected = {
        "SCA1_ATXN1 [CAG_CAA] : valeurs 0–5 hors de toute plage → 'unclassified'",
        "SCA7_ATXN7 [CAG_CAA] : valeurs 20–27 hors de toute plage → 'unclassified'",
        "SCA36_NOP56 [GGCCTG] : valeurs 15–649 hors de toute plage → 'unclassified'",
        "CANVAS_RFC1 [AAGGG] : valeurs 0–199 hors de toute plage → 'unclassified'",
        "CANVAS_RFC1 [AAAGG] : valeurs 0–499 hors de toute plage → 'unclassified'",
        "FXS_FMR1 [CGG] : plages 'premutation' et 'pathogenic_complete' qui se chevauchent sur 200 → 'unclassified'",
        "FXS_FMR1 : clé inconnue 'motif_properties.non_pathogenic_motifs' (ignorée)",
    }
    assert expected <= set(warnings)
    assert len(warnings) == 14


def test_coverage_gaps():
    assert coverage_gaps([(3, 14), (650, None)]) == [(0, 2), (15, 649)]
    assert coverage_gaps([(0, 10), (11, None)]) == []
    assert coverage_gaps([(0, 10)]) == [(11, None)]
    assert coverage_gaps([]) == [(0, None)]


def _mutate(data, fn):
    data = copy.deepcopy(data)
    fn(data)
    return validate_thresholds(data)[0]


def test_unknown_label_is_an_error(thresholds_data):
    def fn(d):
        d["SCA2_ATXN2"]["thresholds"]["CAG_CAA"]["borderline"] = [32, 34]
    assert any("label 'borderline' absent de label_priority" in e for e in _mutate(thresholds_data, fn))


def test_min_greater_than_max_is_an_error(thresholds_data):
    def fn(d):
        d["SCA2_ATXN2"]["thresholds"]["CAG_CAA"]["normal"] = [31, 14]
    assert any("min > max" in e for e in _mutate(thresholds_data, fn))


def test_duplicate_priority_is_an_error(thresholds_data):
    def fn(d):
        d["label_priority"]["ambiguous"] = 0
    assert any("priorité 0 dupliquée" in e for e in _mutate(thresholds_data, fn))


def test_missing_orientation_is_an_error(thresholds_data):
    def fn(d):
        del d["SCA3_ATXN3"]["orientation"]
    assert "SCA3_ATXN3 : 'orientation' manquant" in _mutate(thresholds_data, fn)


def test_invalid_enum_is_an_error(thresholds_data):
    def fn(d):
        d["SCA3_ATXN3"]["repeat_mode"] = "sum"
    assert any("'repeat_mode' invalide" in e for e in _mutate(thresholds_data, fn))


def test_thresholds_on_unknown_group_is_an_error(thresholds_data):
    def fn(d):
        d["SCA3_ATXN3"]["thresholds"]["CAG"] = {"normal": [0, 10]}
    assert any("groupe inconnu 'CAG'" in e for e in _mutate(thresholds_data, fn))


def test_rule_with_unknown_label_is_an_error(thresholds_data):
    def fn(d):
        d["SCA1_ATXN1"]["structure_rules"]["CAG_CAA"][0]["conditions"]["classification"] = "benign"
    assert any("règle 1 : label 'benign'" in e for e in _mutate(thresholds_data, fn))
