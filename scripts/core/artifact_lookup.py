"""
Recherche des fichiers d'un patient dans les archives d'un run (BAM, graphiques).

Les correspondances sont exactes et suivent les conventions de nommage de
tgv_inputs_builder.py : un identifiant patient ne doit jamais sélectionner
les fichiers d'un autre patient (ex. S1 / S10).
"""
import io
import os
import re
import logging
import zipfile

_VCF_SUFFIX_RE = re.compile(r"(\.trgt)?(\.sorted)?\.vcf(\.gz)?$", re.IGNORECASE)


class AmbiguousArtifactError(Exception):
    """Plusieurs fichiers correspondent au même patient : aucun ne doit être choisi."""


def sample_id_from_vcf_name(vcf_name):
    """
    Identifiant patient à partir du nom du VCF.
    Ex : 'S1.trgt.sorted.vcf' → 'S1' ; 'dir/S1.trgt.vcf.gz' → 'S1'.
    """
    return _VCF_SUFFIX_RE.sub("", os.path.basename(vcf_name))


def _zip_names(zip_path):
    with zipfile.ZipFile(zip_path, "r") as z:
        return z.namelist()


def _single(matches, description):
    if len(matches) > 1:
        raise AmbiguousArtifactError(
            f"{len(matches)} fichiers correspondent à {description} : {', '.join(matches)}"
        )
    return matches[0] if matches else None


def _paired_bai(names, bam):
    for candidate in (bam + ".bai", bam[:-4] + ".bai"):
        if candidate in names:
            return candidate
    return None


def _with_bai(zip_path, names, bam):
    if bam is None:
        return None
    bai = _paired_bai(names, bam)
    if bai is None:
        logging.warning(f"BAM '{bam}' found in '{zip_path}' but its index (.bai) is missing.")
        return None
    return zip_path, bam, bai


def find_spanning_bam(zip_path, sample_id):
    """
    BAM spanning d'un patient : '{sample_id}.trgt.spanning.sorted.bam' + index.
    Retourne (zip_path, bam, bai) ou None ; lève AmbiguousArtifactError si doublon.
    """
    names = _zip_names(zip_path)
    target = f"{sample_id}.trgt.spanning.sorted.bam"
    matches = [n for n in names if os.path.basename(n) == target]
    return _with_bai(zip_path, names, _single(matches, f"'{target}'"))


def find_mapped_bam(zip_path, sample_id):
    """
    BAM d'entrée (repeat_reads) d'un patient.

    1. Convention actuelle du builder : '{sample_id}.repeat_reads.bam'.
    2. Anciens runs (noms d'origine des BAM) : l'identifiant doit apparaître
       comme un mot entier dans le nom du fichier (S1 ne correspond pas à S10),
       et un seul BAM doit correspondre.
    """
    names = _zip_names(zip_path)
    bams = [n for n in names if n.endswith(".bam")]

    target = f"{sample_id}.repeat_reads.bam"
    exact = [n for n in bams if os.path.basename(n) == target]
    if exact:
        return _with_bai(zip_path, names, _single(exact, f"'{target}'"))

    token = re.compile(rf"(?<![A-Za-z0-9]){re.escape(sample_id)}(?![A-Za-z0-9])")
    legacy = [n for n in bams if token.search(os.path.basename(n))]
    return _with_bai(zip_path, names, _single(legacy, f"l'identifiant '{sample_id}'"))


def find_plot(zip_path, sample_id, category, trid):
    """
    Graphique TRGT d'un patient pour un locus.

    Archive de run → '{sample_id}_{category}.trvz_alleles.zip' → '{trid}.trvz.svg'.
    Retourne (inner_zip, svg_file) ou None.
    """
    inner_target = f"{sample_id}_{category}.trvz_alleles.zip"
    with zipfile.ZipFile(zip_path, "r") as outer:
        matches = [n for n in outer.namelist() if os.path.basename(n) == inner_target]
        inner_name = _single(matches, f"'{inner_target}'")
        if inner_name is None:
            return None
        with zipfile.ZipFile(io.BytesIO(outer.read(inner_name))) as inner:
            svg = f"{trid}.trvz.svg"
            if svg in inner.namelist():
                return inner_name, svg
    return None


def get_spanning_bam(paths, sample_id):
    """find_spanning_bam sur l'archive 'spanning_bam' du run, si elle existe."""
    zip_path = paths.get("spanning_bam")
    if not zip_path or not os.path.isfile(zip_path):
        return None
    return find_spanning_bam(zip_path, sample_id)


def get_mapped_bam(paths, sample_id):
    """find_mapped_bam sur l'archive 'repeat_reads' du run, si elle existe."""
    zip_path = paths.get("repeat_reads") or paths.get("mapped_bam")
    if not zip_path or not os.path.isfile(zip_path):
        return None
    return find_mapped_bam(zip_path, sample_id)
