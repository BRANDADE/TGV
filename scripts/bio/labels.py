"""
Statuts d'allèle et labels réservés, hors de label_priority (clinical_thresholds.yaml).
"""

# Statut d'appel d'un allèle (scripts/models/allele.py)
CALLED = "called"
NO_CALL = "no_call"   # allèle attendu mais non appelé par TRGT (GT '.')
ABSENT = "absent"     # 2e allèle d'un locus haploïde (ex. chrX, caryotype XY)

# Affichage d'un allèle absent (locus haploïde)
ABSENT_DISPLAY = "-"

# Allèle appelé mais non classable : valeur hors de toute plage du YAML,
# plages qui se chevauchent, ou motifs TRGT absents des groupes cliniques.
# Jamais remplacé par le label de plus basse priorité ("normal").
UNCLASSIFIED = "unclassified"

# Notes expliquant un 'unclassified' (reprises dans les commentaires d'export)
NOTE_OUT_OF_RANGE = "valeur hors des plages du YAML"
NOTE_MOTIF_DISCORDANCE = "motifs TRGT absents des groupes cliniques du YAML"


def note_overlap(labels):
    return f"plages du YAML qui se chevauchent ({' / '.join(labels)})"


def label_score(label, label_priority):
    """
    Score utilisé pour choisir le groupe clinique gagnant.
    'unclassified' se place juste au-dessus du label minimal : un groupe non
    classable n'est jamais masqué par un groupe 'normal', mais tout label
    supérieur (intermédiaire, pathogène…) l'emporte.
    """
    if label == UNCLASSIFIED:
        return min(label_priority.values()) + 0.5
    return label_priority[label]


def is_low_label(label, min_label):
    """Label minimal ou non classable : aucun label clinique significatif."""
    return label in (min_label, UNCLASSIFIED)
