"""F — arrondi des motifs compensés (m) : demis vers le haut, indépendamment de la parité."""
import pytest

from scripts.core.motif_utils import compute_m


@pytest.mark.parametrize("bp, motif_len, expected", [
    (3, 6, 1),    # 0,5 → 1 (round() donnait 0)
    (9, 6, 2),    # 1,5 → 2
    (15, 6, 3),   # 2,5 → 3 (round() donnait 2)
    (1, 3, 0),
    (2, 3, 1),
    (3, 3, 1),
    (0, 3, 0),
    (5, 0, 0),
])
def test_compute_m_rounds_half_up(bp, motif_len, expected):
    assert compute_m(bp, motif_len) == expected
