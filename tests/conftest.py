import os

import pytest
import yaml
from trgt_vcf import vcf_text, write_zip

from scripts.core.analysis import build_run_trids, run_analysis
from scripts.core.vcf_parser import parse_vcf_for_sample

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_thresholds():
    with open(os.path.join(ROOT, "configs", "clinical_thresholds.yaml"), encoding="utf-8") as f:
        return yaml.safe_load(f)


@pytest.fixture
def thresholds_data():
    """Contenu du clinical_thresholds.yaml livré avec TGV (non modifié)."""
    return load_thresholds()


@pytest.fixture
def analyze(tmp_path):
    """
    Analyse un patient synthétique de bout en bout (parsing → classification → affichage).
    Retourne {trid_id: Result}.
    """
    def _analyze(records, thresholds=None, sample="S1", selected=None, aliases=None):
        thresholds = thresholds or load_thresholds()
        vcf_name = f"{sample}.trgt.sorted.vcf"
        zip_path = write_zip(tmp_path / f"{sample}-trgt_vcfs.zip", {vcf_name: vcf_text(records, sample=sample)})

        trids, _, run_trids = build_run_trids(zip_path, thresholds, aliases)
        sample_trids = parse_vcf_for_sample(zip_path, vcf_name, run_trids)
        results, _, _ = run_analysis(
            vcf_name, sample_trids, run_trids, selected or trids,
            thresholds["label_priority"], thresholds.get("low_depth_threshold"),
        )
        return {r.trid: r for r in results}

    return _analyze
