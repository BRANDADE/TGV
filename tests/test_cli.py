"""B3/B4 — CLI autonome (sans interface graphique), idempotente, codes retour exploitables."""
import csv
import os
import subprocess
import sys

import pytest
from trgt_vcf import Call, rc, record, vcf_text, write_zip

import tgv_cli
from scripts.core.marking import mark_pathogenic_genotype

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RFC1_MOTIFS = ["AAAAG", "AAAGG", "AAGGG", "ACAGG", "AGGGC", "AAAGGG"]


def _patient_records(expanded):
    return [
        record("SCA2_ATXN2", ["CTG", "TTG"], [Call(rc("CAG" * 22), depth=20), Call(rc("CAG" * (40 if expanded else 23)))], pos=1000),
        record("SCA17_TBP", ["CAG", "CAA"], [Call("CAG" * 36), Call("CAG" * 38)], pos=2000),
        record("CANVAS_RFC1", RFC1_MOTIFS, [Call("AAAGGG" * 11), Call("AAGGG" * 450)], pos=3000),
        record("HTT", ["CAG"], [Call("CAG" * 17), Call("CAG" * 19)], pos=4000),
    ]


@pytest.fixture
def run_zip(tmp_path):
    members = {
        "S10.trgt.sorted.vcf": vcf_text(_patient_records(True), sample="S10"),
        "S1.trgt.sorted.vcf": vcf_text(_patient_records(False), sample="S1"),
    }
    return write_zip(tmp_path / "RUN42-trgt_vcfs.zip", members)


def _read(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_export_is_idempotent_and_clean(tmp_path, run_zip):
    out = str(tmp_path / "export.tsv")
    assert tgv_cli.main(["--zip", run_zip, "--out", out]) == 0
    first = _read(out)
    assert tgv_cli.main(["--zip", run_zip, "--out", out]) == 0
    second = _read(out)

    # 2 patients × 3 loci du panel Ataxie présents : pas de doublons à la 2e exécution
    assert len(first) == len(second) == 6
    assert {r["sample_id"] for r in second} == {"S1", "S10"}
    assert {r["run_id"] for r in second} == {"RUN42"}

    content = open(out, encoding="utf-8").read()
    assert "⚠" not in content and "●" not in content and "\U0001F534" not in content

    s10_sca2 = next(r for r in second if r["sample_id"] == "S10" and r["Locus"] == "SCA2 (ATXN2)")
    assert s10_sca2["Depth"] == "20 / 30"
    assert s10_sca2["Genotype"] == "22 / 40"
    assert s10_sca2["Classification"] == "normal / pathogenic_complete"
    assert s10_sca2["Comments"] == "Low coverage"


def test_rows_of_other_runs_are_kept(tmp_path, run_zip):
    out = str(tmp_path / "export.tsv")
    assert tgv_cli.main(["--zip", run_zip, "--out", out, "--run-id", "OLD"]) == 0
    assert tgv_cli.main(["--zip", run_zip, "--out", out]) == 0
    rows = _read(out)
    assert len(rows) == 12
    assert {r["run_id"] for r in rows} == {"OLD", "RUN42"}


def test_all_loci_includes_trid_without_underscore(tmp_path, run_zip):
    out = str(tmp_path / "export.tsv")
    assert tgv_cli.main(["--zip", run_zip, "--out", out, "--all-loci"]) == 0
    rows = _read(out)
    assert len(rows) == 8
    htt = next(r for r in rows if r["Locus"] == "HTT")
    assert htt["Classification"] == ""


def test_failing_sample_gives_non_zero_exit_but_writes_others(tmp_path):
    broken = "##fileformat=VCFv4.2\n#CHROM\tPOS\n chr1\t1000\tbroken-line\n"
    zip_path = write_zip(tmp_path / "RUN7-trgt_vcfs.zip", {
        "S1.trgt.sorted.vcf": vcf_text(_patient_records(False), sample="S1"),
        "S2.trgt.sorted.vcf": broken.replace(" chr1", "chr1"),
    })
    out = str(tmp_path / "export.tsv")
    assert tgv_cli.main(["--zip", zip_path, "--out", out]) == tgv_cli.EXIT_FAILURE
    assert {r["sample_id"] for r in _read(out)} == {"S1"}


def test_missing_zip_is_a_configuration_error(tmp_path):
    out = str(tmp_path / "export.tsv")
    assert tgv_cli.main(["--zip", str(tmp_path / "absent.zip"), "--out", out]) == tgv_cli.EXIT_CONFIG
    assert not os.path.exists(out)


def test_cli_runs_without_gui_libraries(tmp_path, run_zip):
    # PySimpleGUI et tkinter rendus inimportables : la CLI doit fonctionner quand même.
    out = str(tmp_path / "export.tsv")
    code = (
        "import sys, runpy\n"
        "sys.modules['PySimpleGUI'] = None\n"
        "sys.modules['tkinter'] = None\n"
        f"sys.argv = ['tgv_cli.py', '--zip', {run_zip!r}, '--out', {out!r}]\n"
        "runpy.run_path('tgv_cli.py', run_name='__main__')\n"
    )
    proc = subprocess.run([sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True, env={**os.environ, "DISPLAY": ""})
    assert proc.returncode == 0, proc.stderr
    assert len(_read(out)) == 6


def test_genotype_marking_uses_whole_motifs():
    marked = mark_pathogenic_genotype("11 (AAAGGG) / 450 (AAGGG)", ["AAGGG"], ui=True)
    assert marked == "11 (AAAGGG) / 450 (●AAGGG)"
