"""A1 — un patient ne doit jamais recevoir les fichiers d'un autre (S1 / S10)."""
import pytest
from trgt_vcf import plots_zip_bytes, write_zip

from scripts.core.artifact_lookup import (
    AmbiguousArtifactError,
    find_mapped_bam,
    find_plot,
    find_spanning_bam,
    get_spanning_bam,
    sample_id_from_vcf_name,
)


@pytest.mark.parametrize("name, expected", [
    ("S1.trgt.sorted.vcf", "S1"),
    ("S1.trgt.vcf", "S1"),
    ("S1.trgt.sorted.vcf.gz", "S1"),
    ("dir/P20.vcf", "P20"),
    ("S1.2024.trgt.sorted.vcf", "S1.2024"),
])
def test_sample_id_from_vcf_name(name, expected):
    assert sample_id_from_vcf_name(name) == expected


def _spanning_zip(tmp_path, samples):
    members = {}
    for s in samples:
        members[f"{s}.trgt.spanning.sorted.bam"] = b"bam-" + s.encode()
        members[f"{s}.trgt.spanning.sorted.bam.bai"] = b"bai-" + s.encode()
    return write_zip(tmp_path / "RUN-spanning_BAM.zip", members)


def test_spanning_bam_exact_match_s1_vs_s10(tmp_path):
    # S10 et S11 sont listés AVANT S1 : l'ancienne recherche par sous-chaîne prenait S10.
    zip_path = _spanning_zip(tmp_path, ["S10", "S11", "S1"])
    _, bam, bai = find_spanning_bam(zip_path, "S1")
    assert bam == "S1.trgt.spanning.sorted.bam"
    assert bai == "S1.trgt.spanning.sorted.bam.bai"


def test_spanning_bam_case_sensitive(tmp_path):
    zip_path = _spanning_zip(tmp_path, ["s1"])
    assert find_spanning_bam(zip_path, "S1") is None


def test_spanning_bam_absent(tmp_path):
    zip_path = _spanning_zip(tmp_path, ["S10"])
    assert find_spanning_bam(zip_path, "S1") is None


def test_spanning_bam_without_index_is_ignored(tmp_path):
    zip_path = write_zip(tmp_path / "z.zip", {
        "S1.trgt.spanning.sorted.bam": b"x",
        "S10.trgt.spanning.sorted.bam.bai": b"x",
    })
    assert find_spanning_bam(zip_path, "S1") is None


def test_spanning_bam_alt_index_name(tmp_path):
    zip_path = write_zip(tmp_path / "z.zip", {
        "S1.trgt.spanning.sorted.bam": b"x",
        "S1.trgt.spanning.sorted.bai": b"x",
    })
    assert find_spanning_bam(zip_path, "S1")[2] == "S1.trgt.spanning.sorted.bai"


def test_spanning_bam_duplicate_raises(tmp_path):
    zip_path = write_zip(tmp_path / "z.zip", {
        "a/S1.trgt.spanning.sorted.bam": b"x",
        "a/S1.trgt.spanning.sorted.bam.bai": b"x",
        "b/S1.trgt.spanning.sorted.bam": b"x",
        "b/S1.trgt.spanning.sorted.bam.bai": b"x",
    })
    with pytest.raises(AmbiguousArtifactError):
        find_spanning_bam(zip_path, "S1")


def test_mapped_bam_current_convention(tmp_path):
    zip_path = write_zip(tmp_path / "z.zip", {
        "S10.repeat_reads.bam": b"x", "S10.repeat_reads.bam.bai": b"x",
        "S1.repeat_reads.bam": b"x", "S1.repeat_reads.bam.bai": b"x",
    })
    assert find_mapped_bam(zip_path, "S1")[1] == "S1.repeat_reads.bam"


def test_mapped_bam_legacy_names_use_word_boundaries(tmp_path):
    zip_path = write_zip(tmp_path / "z.zip", {
        "m84_S10.hifi_reads.mapped.bam": b"x", "m84_S10.hifi_reads.mapped.bam.bai": b"x",
        "m84_S1.hifi_reads.mapped.bam": b"x", "m84_S1.hifi_reads.mapped.bam.bai": b"x",
    })
    assert find_mapped_bam(zip_path, "S1")[1] == "m84_S1.hifi_reads.mapped.bam"
    assert find_mapped_bam(zip_path, "S10")[1] == "m84_S10.hifi_reads.mapped.bam"


def test_mapped_bam_legacy_ambiguous_raises(tmp_path):
    zip_path = write_zip(tmp_path / "z.zip", {
        "run1_S1.mapped.bam": b"x", "run1_S1.mapped.bam.bai": b"x",
        "run2_S1.mapped.bam": b"x", "run2_S1.mapped.bam.bai": b"x",
    })
    with pytest.raises(AmbiguousArtifactError):
        find_mapped_bam(zip_path, "S1")


def test_plot_exact_match_does_not_fall_back_to_other_patient(tmp_path):
    # Seul S10 possède le SVG du locus : S1 ne doit PAS recevoir celui de S10.
    zip_path = write_zip(tmp_path / "RUN-trgt_motifs_allele.zip", {
        "S10_motifs_allele.trvz_alleles.zip": plots_zip_bytes(["SCA1_ATXN1"]),
        "S1_motifs_allele.trvz_alleles.zip": plots_zip_bytes(["SCA2_ATXN2"]),
    })
    assert find_plot(zip_path, "S1", "motifs_allele", "SCA1_ATXN1") is None
    assert find_plot(zip_path, "S1", "motifs_allele", "SCA2_ATXN2") == (
        "S1_motifs_allele.trvz_alleles.zip", "SCA2_ATXN2.trvz.svg"
    )


def test_get_spanning_bam_without_archive():
    assert get_spanning_bam({"spanning_bam": None}, "S1") is None
