"""
Commentaires automatiques associés à un résultat (exports HTML et CLI).
"""
from scripts.bio.labels import NO_CALL

LOW_COVERAGE = "Couverture faible"


def is_low_coverage(result, low_depth_threshold):
    """Vrai si un allèle appelé a une profondeur (SD) inférieure au seuil."""
    if low_depth_threshold is None:
        return False
    try:
        threshold = int(low_depth_threshold)
    except (TypeError, ValueError):
        return False
    for status, depth in ((result.status1, result.depth1_raw), (result.status2, result.depth2_raw)):
        if status != "called":
            continue
        try:
            if int(depth) < threshold:
                return True
        except (TypeError, ValueError):
            continue
    return False


def classification_notes(result):
    """Explications des allèles non classés / non appelés."""
    notes = []
    for idx, note, status in (
        (1, result.classification1_note, result.status1),
        (2, result.classification2_note, result.status2),
    ):
        if note:
            notes.append(f"Allèle {idx} non classé : {note}")
        elif status == NO_CALL:
            notes.append(f"Allèle {idx} non appelé par TRGT")
    return notes


def _pair(v1, v2):
    return f"{'None' if v1 in (None, '') else v1} / {'None' if v2 in (None, '') else v2}"


def override_notes(result):
    """
    Traces des surcharges manuelles (interface) : valeur automatique d'origine.
    Copiées vers le SIL avec la ligne, sans changer le format des colonnes.
    """
    notes = []
    if result.classification1_bio or result.classification2_bio:
        auto = _pair(result.classification1_raw, result.classification2_raw)
        applied = _pair(result.classification1_bio, result.classification2_bio)
        if applied != auto:
            notes.append(f"Classification modifiée manuellement (auto : {auto})")
    if result.genotype1_bio is not None or result.genotype2_bio is not None:
        auto = _pair(result.genotype1_raw, result.genotype2_raw)
        applied = _pair(result.genotype1_bio, result.genotype2_bio)
        if applied != auto:
            notes.append(f"Génotype modifié manuellement (auto : {auto})")
    return notes


def result_comments(result, low_depth_threshold=None):
    """Liste des commentaires automatiques d'un résultat."""
    comments = []
    if is_low_coverage(result, low_depth_threshold):
        comments.append(LOW_COVERAGE)
    comments.extend(classification_notes(result))
    comments.extend(override_notes(result))
    return comments
