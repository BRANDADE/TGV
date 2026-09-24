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


def result_comments(result, low_depth_threshold=None):
    """Liste des commentaires automatiques d'un résultat."""
    comments = []
    if is_low_coverage(result, low_depth_threshold):
        comments.append(LOW_COVERAGE)
    comments.extend(classification_notes(result))
    return comments
