"""C1–C4 — tgv_inputs_builder : paramètres, caryotype, fail reads, reprise par manifeste."""
import argparse
import json
import logging
import os
import zipfile

import pytest

import tgv_inputs_builder as builder


@pytest.fixture
def args(tmp_path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n")
    bed = tmp_path / "panel.bed"
    bed.write_text("chr1\t10\t20\tID=SCA1_ATXN1;MOTIFS=CTG\n")
    return argparse.Namespace(
        trgt="trgt", samtools="samtools", bcftools="bcftools",
        reference=str(ref), bed=str(bed), keep_temp=True, params=builder.DEFAULT_PARAMS_FILE,
    )


def _bam(tmp_path, name):
    path = tmp_path / name
    path.write_bytes(b"BAM")
    return str(path)


def test_default_params_file_is_resolved_from_script_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # répertoire courant sans configs/
    genotype, plot, *_ = builder.load_trgt_params(builder.DEFAULT_PARAMS_FILE)
    assert genotype["preset"] == "targeted"
    assert plot["plot_mode"] == "all"


def test_missing_params_file_stops_the_builder(tmp_path):
    with pytest.raises(SystemExit):
        builder.load_trgt_params(str(tmp_path / "absent.json5"))


def test_missing_karyotype_warns_and_is_traced(tmp_path, caplog):
    bam1, bam2 = _bam(tmp_path, "a.bam"), _bam(tmp_path, "b.bam")
    lst = tmp_path / "samples.tsv"
    lst.write_text(f"S1\t{bam1}\nS2\t{bam2}\txy\n")
    with caplog.at_level(logging.WARNING):
        samples = builder.parse_list_samples(str(lst))
    assert samples["S1"]["karyotype"] == "XX" and samples["S1"]["karyotype_source"] == "default"
    assert samples["S2"]["karyotype"] == "XY" and samples["S2"]["karyotype_source"] == "list"
    assert "No karyotype given for sample 'S1'" in caplog.text


def test_karyotype_file_path_keeps_its_case(tmp_path):
    bam = _bam(tmp_path, "a.bam")
    kfile = tmp_path / "Karyotypes.txt"
    kfile.write_text("S1 XY\n")
    lst = tmp_path / "samples.tsv"
    lst.write_text(f"S1\t{bam}\t{kfile}\n")
    assert builder.parse_list_samples(str(lst))["S1"]["karyotype"] == str(kfile)


def test_fail_reads_column(tmp_path, args):
    bam, fail = _bam(tmp_path, "a.bam"), _bam(tmp_path, "a.fail.bam")
    lst = tmp_path / "samples.tsv"
    lst.write_text(f"S1\t{bam}\tXY\t{fail}\n")
    info = builder.parse_list_samples(str(lst))["S1"]
    assert info["fail_reads"] == fail

    params = {"preset": "targeted", "min-read-quality": 0.98}
    cmd, _ = builder.build_trgt_command("S1", bam, "XY", 4, args, params, str(tmp_path / "out"), fail_reads=fail)
    assert cmd[cmd.index("--fail-reads") + 1] == fail
    assert cmd[cmd.index("--min-read-quality") + 1] == "0.98"
    assert cmd[cmd.index("--preset") + 1] == "targeted"

    cmd, _ = builder.build_trgt_command("S1", bam, "XY", 4, args, params, str(tmp_path / "out"))
    assert "--fail-reads" not in cmd


def test_fail_reads_require_trgt_5_1_and_explicit_quality():
    samples = {"S1": {"fail_reads": "f.bam"}}
    with pytest.raises(SystemExit):
        builder.check_fail_reads_requirements(samples, {"min-read-quality": 0.98}, "5.0.0")
    with pytest.raises(SystemExit):
        builder.check_fail_reads_requirements(samples, {}, "5.1.0")
    builder.check_fail_reads_requirements(samples, {"min-read-quality": 0.98}, "5.1.0")
    builder.check_fail_reads_requirements({"S1": {"fail_reads": None}}, {}, "")


def _fake_genotyping_outputs(output_root, sample_id):
    out = os.path.join(output_root, sample_id)
    os.makedirs(out, exist_ok=True)
    for name in (f"{sample_id}.trgt.sorted.vcf.gz", f"{sample_id}.trgt.sorted.vcf.gz.tbi",
                 f"{sample_id}.trgt.spanning.sorted.bam", f"{sample_id}.trgt.spanning.sorted.bam.bai"):
        with open(os.path.join(out, name), "wb") as f:
            f.write(b"x")


def test_resume_only_with_identical_manifest(tmp_path, args, monkeypatch):
    output_root = str(tmp_path / "inputs_RUN")
    bam = _bam(tmp_path, "a.bam")
    info = {"bam_path": bam, "karyotype": "XX", "fail_reads": None}
    params = {"preset": "targeted"}
    fingerprint = builder.genotype_fingerprint("S1", info, args, params, "5.1.0", "abc")

    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        raise builder.subprocess.CalledProcessError(1, cmd, stderr="stop")

    monkeypatch.setattr(builder.subprocess, "run", fake_run)
    _fake_genotyping_outputs(output_root, "S1")
    open(bam + ".bai", "wb").close()

    # Sorties existantes sans manifeste → nouveau génotypage
    with pytest.raises(builder.subprocess.CalledProcessError):
        builder.run_trgt_genotype("S1", bam, "XX", 1, args, params, output_root, fingerprint=fingerprint)
    assert calls and calls[-1][1] == "genotype"

    # Manifeste identique → réutilisation, aucun appel à TRGT
    builder.write_manifest(builder.sample_manifest_path(output_root, "S1"), fingerprint)
    calls.clear()
    builder.run_trgt_genotype("S1", bam, "XX", 1, args, params, output_root, fingerprint=fingerprint)
    assert calls == []

    # BED modifié → nouveau génotypage
    changed = builder.genotype_fingerprint("S1", info, args, params, "5.1.0", "def")
    with pytest.raises(builder.subprocess.CalledProcessError):
        builder.run_trgt_genotype("S1", bam, "XX", 1, args, params, output_root, fingerprint=changed)
    assert calls


def test_fingerprint_ignores_threads_and_verbosity(args):
    info = {"bam_path": "a.bam", "karyotype": "XX", "fail_reads": None}
    a = builder.genotype_fingerprint("S1", info, args, {"preset": "targeted", "verbose": 0}, "5.1.0", "x")
    b = builder.genotype_fingerprint("S1", info, args, {"preset": "targeted", "verbose": 2}, "5.1.0", "x")
    c = builder.genotype_fingerprint("S1", info, args, {"preset": "wgs"}, "5.1.0", "x")
    assert a == b and a != c


def test_aggregation_names_input_bams_per_sample_and_adds_manifest(tmp_path):
    output_root = str(tmp_path / "inputs_RUN")
    os.makedirs(os.path.join(output_root, "S1"))
    bam = _bam(tmp_path, "m84_run_S1.hifi_reads.bam")
    open(bam + ".bai", "wb").close()
    import gzip
    with gzip.open(os.path.join(output_root, "S1", "S1.trgt.sorted.vcf.gz"), "wt") as f:
        f.write("##fileformat=VCFv4.2\n")

    manifest = {"run_name": "RUN", "samples": {"S1": {"karyotype": "XY"}}}
    builder.create_global_aggregations(output_root, "RUN", {"S1": {"bam_path": bam}}, skip_plots=True,
                                       run_manifest=manifest)

    with zipfile.ZipFile(os.path.join(output_root, "RUN-repeat_reads.zip")) as z:
        assert sorted(z.namelist()) == ["S1.repeat_reads.bam", "S1.repeat_reads.bam.bai"]
    with zipfile.ZipFile(os.path.join(output_root, "RUN-trgt_vcfs.zip")) as z:
        assert json.loads(z.read("run_manifest.json"))["samples"]["S1"]["karyotype"] == "XY"
    assert os.path.exists(os.path.join(output_root, "RUN-manifest.json"))


# --- Rapport QC : métriques de lectures et profondeur par allèle -----------------

def _sam(flag, seq, qual):
    return "\t".join(["r", str(flag), "chr1", "1", "60", f"{len(seq)}M", "*", "0", "0", seq, qual]) + "\n"


def test_read_qv_averages_error_probabilities():
    # Q40 partout → Q40 ; un mélange Q10/Q40 est dominé par les erreurs (≈ Q13), pas Q25
    assert round(builder.read_qv("I" * 10)) == 40
    assert round(builder.read_qv("+" * 5 + "I" * 5)) == 13
    assert builder.read_qv("*") is None


def test_hist_median_uses_sorted_values():
    assert builder.hist_median({5000: 1, 4800: 1, 5200: 1}) == 5000
    assert builder.hist_median({}) == 0


def test_summarise_reads_counts_all_primary_reads():
    lines = [_sam(0, "A" * 100, "I" * 100), _sam(16, "C" * 300, "I" * 300),
             _sam(0x100, "G" * 50, "I" * 50),          # secondaire : ignorée
             _sam(0x800, "T" * 70, "I" * 70)]          # supplémentaire : ignorée
    m = builder.summarise_reads(lines * 6000)            # > 10 000 lectures : aucun échantillonnage
    assert m["num_reads"] == 12000
    assert m["total_bases"] == 6000 * 400
    assert m["median_read_length"] == 100
    assert m["median_read_qv"] == 40


def test_qc_depth_sums_homozygous_calls():
    assert builder.qc_depth([43, 43], "1/1") == 86
    assert builder.qc_depth([100, 8], "0/1") == 8
    assert builder.qc_depth([40], "1") == 40
    assert builder.qc_depth([], ".") == 0


def test_parse_vcf_for_qc_applies_tgv_rule(tmp_path):
    import gzip
    vcf = tmp_path / "s.trgt.sorted.vcf.gz"
    rows = [
        ("FXS_FMR1", "1/1", "43,43"),     # homozygote : 86 lectures → non signalé
        ("FRDA_FXN", "0/1", "100,8"),     # allèle expansé à 8 lectures → signalé
    ]
    with gzip.open(vcf, "wt") as f:
        f.write("##fileformat=VCFv4.2\n")
        for trid, gt, sd in rows:
            f.write("\t".join(["chr1", "1", ".", "A", "C", ".", ".", f"TRID={trid}",
                               "GT:SD", f"{gt}:{sd}"]) + "\n")
    details, summary = builder.parse_vcf_for_qc(str(vcf), dp_threshold=50)
    assert details["FXS_FMR1"]["min_dp"] == 86
    assert [f["repeat_id"] for f in summary["flagged_loci"]] == ["FRDA_FXN"]


def test_default_params_file_declares_qc_thresholds():
    *_, qc = builder.load_trgt_params(builder.DEFAULT_PARAMS_FILE)
    assert qc == {"min_allele_depth": 50, "min_mean_coverage": 10, "min_pct_pass_loci": 90}
