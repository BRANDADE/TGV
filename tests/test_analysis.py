"""B5 — l'analyse ne plante pas quand des loci sélectionnés sont absents du VCF du patient."""
from trgt_vcf import Call, vcf_text, record, write_zip

from scripts.core.analysis import build_run_trids, run_analysis
from scripts.core.vcf_parser import parse_vcf_for_sample


def test_selected_loci_absent_from_patient_are_reported(tmp_path, thresholds_data):
    # Le 1er VCF (référence du catalogue) contient SCA3 ; celui de S2 non.
    zip_path = write_zip(tmp_path / "RUN-trgt_vcfs.zip", {
        "S1.trgt.sorted.vcf": vcf_text([
            record("SCA17_TBP", ["CAG", "CAA"], [Call("CAG" * 36), Call("CAG" * 38)], pos=1000),
            record("SCA3_ATXN3", ["CTG", "TTG"], [Call("CTG" * 20), Call("CTG" * 22)], pos=2000),
        ]),
        "S2.trgt.sorted.vcf": vcf_text([
            record("SCA17_TBP", ["CAG", "CAA"], [Call("CAG" * 36), Call("CAG" * 38)], pos=1000),
        ], sample="S2"),
    })
    trids, _, run_trids = build_run_trids(zip_path, thresholds_data)
    sample_trids = parse_vcf_for_sample(zip_path, "S2.trgt.sorted.vcf", run_trids)

    results, discordances, missing = run_analysis(
        "S2.trgt.sorted.vcf", sample_trids, run_trids, ["SCA3_ATXN3", "SCA17_TBP", "UNKNOWN"],
        thresholds_data["label_priority"], thresholds_data["low_depth_threshold"],
    )
    assert [r.trid for r in results] == ["SCA17_TBP"]
    assert missing == ["SCA3_ATXN3", "UNKNOWN"]
    assert discordances == []


def test_repeated_analysis_is_stable(tmp_path, thresholds_data):
    zip_path = write_zip(tmp_path / "RUN-trgt_vcfs.zip", {
        "S1.trgt.sorted.vcf": vcf_text([
            record("SCA17_TBP", ["CAG", "CAA"], [Call("CAG" * 52), Call("CAG" * 36)]),
        ]),
    })
    _, _, run_trids = build_run_trids(zip_path, thresholds_data)
    sample_trids = parse_vcf_for_sample(zip_path, "S1.trgt.sorted.vcf", run_trids)
    args = (["SCA17_TBP"], thresholds_data["label_priority"], thresholds_data["low_depth_threshold"])

    first, _, _ = run_analysis("S1.trgt.sorted.vcf", sample_trids, run_trids, *args)
    second, _, _ = run_analysis("S1.trgt.sorted.vcf", sample_trids, run_trids, *args)
    assert first[0] is second[0]
    assert second[0].display_export.genotype == "36 / 52"
