import os

import pytest
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture
def thresholds_data():
    """Contenu du clinical_thresholds.yaml livré avec TGV (non modifié)."""
    with open(os.path.join(ROOT, "configs", "clinical_thresholds.yaml"), encoding="utf-8") as f:
        return yaml.safe_load(f)
