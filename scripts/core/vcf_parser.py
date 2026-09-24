import zipfile
import logging

from scripts.models.sample import Sample
from scripts.models.allele import Allele, AlleleSequence
from scripts.bio.labels import CALLED, NO_CALL, ABSENT

from scripts.core.sequence_utils import reverse_complement, rc_segmentation, convert_mc, convert_ms
from scripts.core.segmentation_interruptions import find_interruptions, extract_interruption_sequences, segmentation_complete


def parse_info(info_str):
    """Parse the INFO field of a VCF line."""
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
    return info


def split_values(value):
    """
    Sépare un champ FORMAT par allèle ('x,y' → ['x', 'y']).
    Les valeurs manquantes VCF '.' deviennent des chaînes vides.
    """
    if value is None:
        return []
    return ["" if v == "." else v for v in value.split(",")]


def parse_gt(gt):
    """
    Indices VCF des allèles appelés, dans l'ordre des champs FORMAT.

    TRGT écrit les champs FORMAT (AL, ALLR, SD, MC, MS, AP, AM) dans son ordre
    interne des allèles, allèle de référence en premier : les indices VCF sont
    donc croissants dans cet ordre. Un GT phasé (ex. '1|0') n'est qu'une
    permutation qui ne réordonne PAS les champs FORMAT.
    Réf. : TRGT docs/vcf_files.md ; src/trgt/writers/write_vcf.rs (set_gt) ;
           src/trgt/workflows/tr.rs (reference allele first).

    Retourne (indices triés, ploïdie). Ex : '1|0' → ([0, 1], 2) ; '1' → ([1], 1) ;
    '.' → ([], 0) (no-call).
    """
    if not gt or gt == ".":
        return [], 0

    parts = gt.replace("|", "/").split("/")
    indices = []
    for part in parts:
        try:
            indices.append(int(part))
        except ValueError:
            continue  # allèle manquant '.'
    return sorted(indices), len(parts)


def allele_sequence(index, ref, alt_list):
    """Séquence de l'allèle d'indice VCF donné, sans la base de padding TRGT."""
    if index == 0:
        seq = ref
    elif index - 1 < len(alt_list):
        seq = alt_list[index - 1]
    else:
        logging.warning(
            f"VCF format anomaly: genotype index '{index}' points to a non-existent "
            f"alternative allele (ALT list len={len(alt_list)}). Defaulting to REF."
        )
        seq = ref
    # TRGT : toujours une base de padding en tête → on l'enlève
    return seq[1:] if len(seq) > 1 else ""


def _value(values, i):
    return values[i] if i < len(values) else ""


def empty_allele(status):
    """Allèle non appelé (no_call) ou absent (locus haploïde)."""
    allele = Allele(
        size="", size_range="", depth="", purity="", methylation="",
        sequence=AlleleSequence("", "", ""),
        status=status,
    )
    allele.sequence.interruptions = []
    allele.sequence.segmentation_complete = ""
    return allele


def build_called_allele(cons, fields, i, motifs, orientation, motif_names=None):
    """
    Construit un allèle appelé à partir de sa séquence et de ses champs FORMAT (position i).

    motifs      : motifs du catalogue, orientés (servent à la segmentation RC) ;
    motif_names : noms affichés/classés (cadre de lecture du YAML), même ordre.
    """
    motif_names = motif_names or motifs
    mc_raw = _value(fields["MC"], i)
    ms_raw = _value(fields["MS"], i)

    # ORIENTATION RC → recalcul de la séquence et de la segmentation
    if orientation == "rc":
        cons = reverse_complement(cons)
        ms_raw = rc_segmentation(cons, ",".join(motifs))

    # Les comptes (MC) et les coordonnées (MS) restent ceux de TRGT ; seuls les noms
    # des motifs sont ramenés au cadre de lecture du YAML (ex. GCA → CAG).
    mc = convert_mc(motif_names, mc_raw)
    ms = convert_ms(motif_names, ms_raw)

    allele = Allele(
        size=_value(fields["AL"], i),
        size_range=_value(fields["ALLR"], i),
        depth=_value(fields["SD"], i),
        purity=_value(fields["AP"], i),
        methylation=_value(fields["AM"], i),
        sequence=AlleleSequence(cons, mc, ms),
        status=CALLED,
    )

    intervals = find_interruptions(ms)
    inter_seqs = extract_interruption_sequences(cons, intervals)
    allele.sequence.interruptions = inter_seqs
    allele.sequence.segmentation_complete = segmentation_complete(ms, intervals, inter_seqs)
    return allele


def parse_vcf_for_sample(zip_path, vcf_filename, global_trids):
    """
    global_trids : dict { trid_id: TRID } venant de run.trids
    Retourne : dict { trid_id: Sample }
    """
    samples = {}

    logging.info(f"Opening VCF file '{vcf_filename}' inside archive '{zip_path}'")

    with zipfile.ZipFile(zip_path, "r") as z:
        # Sécurité : vérifier que le fichier VCF existe bien dans l'archive
        if vcf_filename not in z.namelist():
            logging.error(f"VCF file '{vcf_filename}' was not found inside TRGT ZIP archive '{zip_path}'.")
            raise FileNotFoundError(f"VCF file '{vcf_filename}' not found in ZIP.")

        with z.open(vcf_filename) as f:
            for raw in f:
                line = raw.decode("utf-8").strip()
                if not line or line.startswith("#"):
                    continue

                cols = line.split("\t")
                chrom, pos, vid, ref, alt, qual, flt, info_str, fmt, sample_data = cols

                info = parse_info(info_str)
                trid_id = info.get("TRID")
                if not trid_id:
                    continue

                # On ne traite que les TRIDs connus globalement
                trid_global = global_trids.get(trid_id)
                if not trid_global:
                    continue

                # -----------------------------
                # ORIENTATION CLINIQUE (FW / RC)
                # -----------------------------
                orientation = "fw"
                motifs = list(trid_global.motifs)  # copie

                if trid_global.clinical and trid_global.clinical.orientation:
                    if trid_global.clinical.orientation.lower() == "rc":
                        orientation = "rc"
                        motifs = [reverse_complement(m) for m in motifs]

                motif_frame = getattr(trid_global, "motif_frame", None) or {}
                motif_names = [motif_frame.get(m, m) for m in motifs]

                # -----------------------------
                # FORMAT → extraction TRGT (valeurs dans l'ordre des allèles TRGT)
                # -----------------------------
                data = dict(zip(fmt.split(":"), sample_data.split(":")))
                fields = {key: split_values(data.get(key)) for key in ("AL", "ALLR", "SD", "MC", "MS", "AP", "AM")}

                sample_obj = Sample(name=vcf_filename)
                sample_obj.gt_raw = data.get("GT")

                indices, ploidy = parse_gt(sample_obj.gt_raw)
                alt_list = [] if alt == "." else alt.split(",")

                alleles = []
                for i, index in enumerate(indices[:2]):
                    cons = allele_sequence(index, ref, alt_list)
                    alleles.append(build_called_allele(cons, fields, i, motifs, orientation, motif_names))

                # Complément à deux allèles : haploïde → absent ; sinon non appelé
                while len(alleles) < 2:
                    status = ABSENT if (ploidy == 1 and len(alleles) == 1) else NO_CALL
                    alleles.append(empty_allele(status))

                sample_obj.allele1, sample_obj.allele2 = alleles
                samples[trid_id] = sample_obj

    # Bilan de fin de parsing pour l'audit de qualité
    if not samples:
        logging.warning(f"No genomic loci were successfully parsed from VCF file '{vcf_filename}'.")
    else:
        logging.info(f"Successfully parsed {len(samples)} loci from patient VCF '{vcf_filename}'.")

    return samples
