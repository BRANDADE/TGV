"""
Statuts d'allèle et labels réservés, hors de label_priority (clinical_thresholds.yaml).
"""

# Statut d'appel d'un allèle (scripts/models/allele.py)
CALLED = "called"
NO_CALL = "no_call"   # allèle attendu mais non appelé par TRGT (GT '.')
ABSENT = "absent"     # 2e allèle d'un locus haploïde (ex. chrX, caryotype XY)

# Affichage d'un allèle absent (locus haploïde)
ABSENT_DISPLAY = "-"
