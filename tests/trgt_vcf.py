"""
Fabrique de VCF synthétiques au format TRGT 5.x (tests uniquement).

Le format reproduit src/trgt/writers/write_vcf.rs de TRGT :
- champs FORMAT GT:AL:ALLR:SD:MC:MS:AP:AM, dans l'ordre des allèles TRGT ;
- allèle de référence en premier (src/trgt/workflows/tr.rs), donc GT non phasé croissant ;
- un GT phasé (ex. "1|0") est une permutation qui ne réordonne PAS les champs FORMAT ;
- base de padding en tête de REF/ALT ;
- no-call : GT "." et tous les champs FORMAT à ".".
"""
import io
import zipfile
from dataclasses import dataclass
from typing import Optional

FORMAT = "GT:AL:ALLR:SD:MC:MS:AP:AM"
_COMP = str.maketrans("ACGTN", "TGCAN")


def rc(seq):
    """Reverse-complement."""
    return seq.translate(_COMP)[::-1]


def segment(seq, motifs):
    """Segmentation gloutonne : renvoie (comptes par motif, chaîne MS, bases couvertes)."""
    order = sorted(range(len(motifs)), key=lambda i: len(motifs[i]), reverse=True)
    counts = [0] * len(motifs)
    spans = []
    covered = 0
    i = 0
    while i < len(seq):
        for idx in order:
            motif = motifs[idx]
            if seq.startswith(motif, i):
                start = i
                while seq.startswith(motif, i):
                    i += len(motif)
                    counts[idx] += 1
                spans.append(f"{idx}({start}-{i})")
                covered += i - start
                break
        else:
            i += 1
    return counts, "_".join(spans), covered


@dataclass
class Call:
    """Un allèle appelé par TRGT (séquence du brin +, sans base de padding)."""
    seq: str
    depth: int = 30
    methylation: Optional[float] = None
    mc: Optional[str] = None  # surcharge de MC (sinon calculé)
    ms: Optional[str] = None  # surcharge de MS (sinon calculé ; "." accepté)


def _vcf_digits(calls, ref):
    """Reproduit VcfWriter::set_gt : numérotation VCF des allèles et liste ALT."""
    alts = []
    digits = []
    for call in calls:
        if call.seq == ref:
            digits.append(0)
            continue
        if call.seq not in alts:
            alts.append(call.seq)
        digits.append(alts.index(call.seq) + 1)
    return digits, alts


def record(trid, motifs, calls, *, ref=None, gt=None, chrom="chr1", pos=1000, pad="A", struc="<TR>"):
    """
    Construit une ligne de données VCF TRGT.

    calls : allèles dans l'ordre TRGT (référence en premier) ; 1 allèle = haploïde,
            liste vide = no-call.
    gt    : GT écrit tel quel (ex. "1|0") ; par défaut, GT non phasé dérivé des allèles.
    """
    if ref is None:
        ref = calls[0].seq if calls else ""

    info = f"TRID={trid};END={pos + len(ref)};MOTIFS={','.join(motifs)};STRUC={struc}"

    if not calls:
        sample = ":".join(["."] * len(FORMAT.split(":")))
        cols = [chrom, str(pos), ".", pad + ref, ".", ".", ".", info, FORMAT, sample]
        return "\t".join(cols)

    digits, alts = _vcf_digits(calls, ref)
    if gt is None:
        gt = "/".join(str(d) for d in digits)

    al, allr, sd, mc, ms, ap, am = [], [], [], [], [], [], []
    for call in calls:
        counts, ms_auto, covered = segment(call.seq, motifs)
        length = len(call.seq)
        al.append(str(length))
        allr.append(f"{length}-{length}")
        sd.append(str(call.depth))
        mc.append(call.mc if call.mc is not None else "_".join(str(c) for c in counts))
        ms.append(call.ms if call.ms is not None else ms_auto)
        ap.append(f"{covered / length:.6g}" if length else ".")
        am.append("." if call.methylation is None else f"{call.methylation:.2f}")

    sample = ":".join([
        gt, ",".join(al), ",".join(allr), ",".join(sd),
        ",".join(mc), ",".join(ms), ",".join(ap), ",".join(am),
    ])
    alt_field = ",".join(pad + a for a in alts) if alts else "."
    cols = [chrom, str(pos), ".", pad + ref, alt_field, ".", ".", info, FORMAT, sample]
    return "\t".join(cols)


def vcf_text(records, sample="S1", version="5.1.0", repeats="catalog.bed"):
    """Assemble un VCF complet avec l'en-tête TRGT (##trgtVersion, ##trgtCommand)."""
    header = [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=TRID,Number=1,Type=String,Description="Tandem repeat ID">',
        f"##trgtVersion={version}",
        f"##trgtCommand=trgt genotype --genome ref.fa --repeats {repeats} "
        f"--reads {sample}.bam --output-prefix {sample}.trgt --preset targeted",
        "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]),
    ]
    return "\n".join(header + list(records)) + "\n"


def write_zip(path, members):
    """Écrit un ZIP {nom: contenu (str ou bytes)} et renvoie son chemin (str)."""
    with zipfile.ZipFile(path, "w") as z:
        for name, content in members.items():
            z.writestr(name, content)
    return str(path)


def plots_zip_bytes(trids):
    """Contenu d'un ZIP *.trvz_alleles.zip contenant un SVG par TRID."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for trid in trids:
            z.writestr(f"{trid}.trvz.svg", f"<svg><!-- {trid} --></svg>")
    return buf.getvalue()
