"""Vérifie que la fabrique de VCF produit un format lisible par le parseur TGV."""
from trgt_vcf import Call, record, rc, segment, vcf_text, write_zip

from scripts.core.vcf_parser import parse_vcf_for_sample
from scripts.models.trid import TRID


def test_segment_counts_and_spans():
    counts, ms, covered = segment("CAGCAGCATCAG", ["CAG", "CAT"])
    assert counts == [3, 1]
    assert ms == "0(0-6)_1(6-9)_0(9-12)"
    assert covered == 12


def test_record_format_matches_trgt_layout():
    line = record("X_Y", ["CAG"], [Call("CAG" * 10), Call("CAG" * 20)])
    cols = line.split("\t")
    assert cols[3] == "A" + "CAG" * 10
    assert cols[4] == "A" + "CAG" * 20
    assert cols[8] == "GT:AL:ALLR:SD:MC:MS:AP:AM"
    assert cols[9].startswith("0/1:30,60:30-30,60-60:30,30:10,20:0(0-30),0(0-60):")


def test_no_call_record():
    line = record("X_Y", ["CAG"], [], ref="CAG" * 10)
    assert line.split("\t")[9] == ".:.:.:.:.:.:.:."


def test_parser_reads_synthetic_vcf(tmp_path):
    line = record("SCA2_ATXN2", ["CTG", "TTG"], [Call(rc("CAG" * 22)), Call(rc("CAG" * 40))])
    zip_path = write_zip(tmp_path / "RUN-trgt_vcfs.zip", {"S1.trgt.sorted.vcf": vcf_text([line])})

    trid = TRID("SCA2_ATXN2")
    trid.motifs = ["CTG", "TTG"]
    samples = parse_vcf_for_sample(zip_path, "S1.trgt.sorted.vcf", {"SCA2_ATXN2": trid})

    sample = samples["SCA2_ATXN2"]
    assert sample.allele1.size == "66"
    assert sample.allele2.size == "120"
