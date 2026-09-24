"""
Pipeline d'analyse sans interface graphique, partagé par la CLI et l'interface.
"""
import logging

from scripts.models.trid import TRID
from scripts.models.dto import AnalysisInput

from scripts.core.trid_detector import autodetect_trids, make_readable_name
from scripts.core.sequence_utils import reverse_complement, motif_frame_map
from scripts.core.orchestrator import process_clinical, process_result, process_display

from scripts.bio.clinical_config_builder import build_clinical_config


def clinical_key_for(trid_id, thresholds_data, aliases):
    """Bloc de clinical_thresholds.yaml d'un TRID : le TRID lui-même, sinon son alias."""
    if isinstance(thresholds_data.get(trid_id), dict):
        return trid_id
    alias = (aliases or {}).get(trid_id)
    if alias and isinstance(thresholds_data.get(alias), dict):
        return alias
    return None


def build_run_trids(zip_path, thresholds_data, aliases=None):
    """
    TRID globaux du run : informations statiques (premier VCF) + configuration clinique (YAML).

    aliases : {TRID du catalogue: bloc YAML} (configs/trid_aliases.yaml).
    Retourne (liste triée des TRID, {nom lisible: TRID}, {trid_id: TRID}).
    """
    trids, _, _, static_info = autodetect_trids(zip_path)

    run_trids = {}
    diseases = {}
    for trid_id in trids:
        t = TRID(trid_id)

        info = static_info.get(trid_id)
        if info:
            t.chrom = info["chrom"]
            t.start = info["start"]
            t.end = info["end"]
            t.motifs = info["motifs"]

        t.clinical_key = clinical_key_for(trid_id, thresholds_data, aliases)
        if t.clinical_key:
            if t.clinical_key != trid_id:
                logging.info(f"Locus '{trid_id}' uses clinical thresholds of '{t.clinical_key}' (trid_aliases.yaml).")

            t.clinical = build_clinical_config(t.clinical_key, thresholds_data)
            orientation = (t.clinical.orientation or "").lower()
            oriented = list(t.motifs)
            if orientation == "rc":
                t.motifs_rc = [reverse_complement(m) for m in t.motifs]
                oriented = list(t.motifs_rc)

            clinical_motifs = {m for g in t.clinical.groups.values() for m in g.motifs}
            clinical_motifs |= set(t.clinical.pathogenic_motifs)
            t.motif_frame = motif_frame_map(oriented, clinical_motifs)
            if t.motif_frame:
                logging.warning(
                    f"Locus '{trid_id}': catalog motif(s) renamed to the clinical reading frame "
                    f"{t.motif_frame} (same repeat unit); TRGT motif counts are kept."
                )

        run_trids[trid_id] = t
        diseases[make_readable_name(t.clinical_key or trid_id)] = trid_id

    return trids, diseases, run_trids


def resolve_panel(requested, run_trids):
    """
    TRID du run correspondant à un panel (identifiants YAML ou TRID du catalogue).
    Retourne (TRID présents dans l'ordre du panel, entrées du panel introuvables).
    """
    by_key = {}
    for trid_id, t in run_trids.items():
        by_key.setdefault(trid_id, trid_id)
        if t.clinical_key:
            by_key.setdefault(t.clinical_key, trid_id)

    found, missing = [], []
    for name in requested:
        trid_id = by_key.get(name)
        if trid_id is None:
            missing.append(name)
        elif trid_id not in found:
            found.append(trid_id)
    return found, missing


def run_analysis(sample_name, sample_trids, run_trids, selected_trids,
                 label_priority, low_depth_threshold, paths=None):
    """
    Analyse clinique d'un patient pour les TRID sélectionnés.

    sample_trids : {trid_id: Sample} issu de parse_vcf_for_sample.

    Retourne (results, discordances, missing) :
      - results      : liste de Result, dans l'ordre de selected_trids ;
      - discordances : TRID dont les motifs TRGT ne correspondent pas au YAML ;
      - missing      : TRID sélectionnés absents du VCF du patient.
    """
    present = [t for t in selected_trids if t in sample_trids and t in run_trids]
    missing = [t for t in selected_trids if t not in present]

    analysis_input = AnalysisInput(
        sample_name=sample_name,
        trids={t: run_trids[t] for t in present},
        samples={t: sample_trids[t] for t in present},
        ordered_trids=present,
        paths=paths or {},
        label_priority=label_priority,
    )

    discordances = process_clinical(analysis_input)
    process_result(analysis_input)

    results = []
    for trid_id in present:
        result = sample_trids[trid_id].result
        process_display(result, run_trids[trid_id].clinical, low_depth_threshold)
        results.append(result)

    return results, discordances, missing
