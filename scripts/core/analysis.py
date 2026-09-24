"""
Pipeline d'analyse sans interface graphique, partagé par la CLI et l'interface.
"""
from scripts.models.trid import TRID
from scripts.models.dto import AnalysisInput

from scripts.core.trid_detector import autodetect_trids
from scripts.core.sequence_utils import reverse_complement
from scripts.core.orchestrator import process_clinical, process_result, process_display

from scripts.bio.clinical_config_builder import build_clinical_config


def build_run_trids(zip_path, thresholds_data):
    """
    TRID globaux du run : informations statiques (premier VCF) + configuration clinique (YAML).
    Retourne (liste triée des TRID, {nom lisible: TRID}, {trid_id: TRID}).
    """
    trids, _, diseases, static_info = autodetect_trids(zip_path)

    run_trids = {}
    for trid_id in trids:
        t = TRID(trid_id)

        info = static_info.get(trid_id)
        if info:
            t.chrom = info["chrom"]
            t.start = info["start"]
            t.end = info["end"]
            t.motifs = info["motifs"]

        if isinstance(thresholds_data.get(trid_id), dict):
            t.clinical = build_clinical_config(trid_id, thresholds_data)
            orientation = (t.clinical.orientation or "").lower()
            if orientation == "rc":
                t.motifs_rc = [reverse_complement(m) for m in t.motifs]

        run_trids[trid_id] = t

    return trids, diseases, run_trids


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
