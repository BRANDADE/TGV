"""
Validation de clinical_thresholds.yaml au chargement.

- Erreurs : configuration inutilisable ou ambiguë → le chargement est refusé.
- Avertissements : trous, chevauchements, clés inconnues → le chargement continue ;
  à l'exécution, une valeur concernée est classée 'unclassified' (jamais 'normal'
  par défaut). Voir scripts/bio/clinical_classifier.py.
"""
import re

GLOBAL_KEYS = {"label_priority", "low_depth_threshold"}

LOCUS_KEYS = {
    "source", "classification_mode", "repeat_mode", "genotype_display",
    "orientation", "motif_properties", "thresholds", "structure_rules",
}

MOTIF_PROPERTY_KEYS = {
    "pathogenic_motifs", "protective_motifs", "uncertain_motifs", "pure_only", "motif_groups",
}

ENUMS = {
    "classification_mode": {"simple", "structural"},
    "repeat_mode": {"sum_with_interruptions", "sum_without_interruptions"},
    "genotype_display": {"pathogenic_only", "pathogenic_with_motif", "full_with_others"},
    "orientation": {"fw", "rc"},  # comparaison insensible à la casse
}

_MOTIF_RE = re.compile(r"^[ACGTN]+$")


def _is_int(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _check_range(bounds, where, errors):
    """[min, max] entiers ou null, min <= max. Retourne (min, max) ou None si invalide."""
    if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
        errors.append(f"{where} : plage invalide {bounds!r} (attendu [min, max])")
        return None
    lo, hi = bounds
    for v in (lo, hi):
        if v is not None and not _is_int(v):
            errors.append(f"{where} : borne non entière {v!r}")
            return None
    if lo is not None and hi is not None and lo > hi:
        errors.append(f"{where} : min > max ({lo} > {hi})")
        return None
    return lo, hi


def _format_gap(lo, hi):
    if hi is None:
        return f"≥ {lo}"
    return f"{lo}" if lo == hi else f"{lo}–{hi}"


def coverage_gaps(ranges):
    """Intervalles entiers de [0, +∞[ non couverts par les plages (bornes inclusives)."""
    gaps = []
    cursor = 0  # plus petite valeur non encore couverte
    for lo, hi in sorted(ranges, key=lambda r: (r[0] if r[0] is not None else 0)):
        lo = lo if lo is not None else 0
        if lo > cursor:
            gaps.append((cursor, lo - 1))
        if hi is None:
            return gaps
        cursor = max(cursor, hi + 1)
    gaps.append((cursor, None))
    return gaps


def overlaps(labelled_ranges):
    """Paires de labels différents dont les plages se chevauchent : [(l1, l2, lo, hi)]."""
    found = []
    items = list(labelled_ranges.items())
    for i, (l1, (a_lo, a_hi)) in enumerate(items):
        for l2, (b_lo, b_hi) in items[i + 1:]:
            lo = max(a_lo if a_lo is not None else 0, b_lo if b_lo is not None else 0)
            hi_candidates = [h for h in (a_hi, b_hi) if h is not None]
            hi = min(hi_candidates) if hi_candidates else None
            if hi is None or lo <= hi:
                found.append((l1, l2, lo, hi))
    return found


def _validate_label_priority(data, errors):
    label_priority = data.get("label_priority")
    if not isinstance(label_priority, dict) or not label_priority:
        errors.append("label_priority : section absente ou vide")
        return {}

    by_priority = {}
    for label, prio in label_priority.items():
        if not _is_int(prio):
            errors.append(f"label_priority : priorité non entière pour '{label}' ({prio!r})")
            continue
        by_priority.setdefault(prio, []).append(label)
    for prio, labels in sorted(by_priority.items()):
        if len(labels) > 1:
            errors.append(f"label_priority : priorité {prio} dupliquée ({', '.join(labels)})")
    return label_priority


def _validate_locus(trid, block, label_priority, errors, warnings):
    for key in sorted(set(block) - LOCUS_KEYS):
        warnings.append(f"{trid} : clé inconnue '{key}' (ignorée)")

    for key, allowed in ENUMS.items():
        value = block.get(key)
        normalized = value.lower() if (key == "orientation" and isinstance(value, str)) else value
        if value is None:
            errors.append(f"{trid} : '{key}' manquant")
        elif normalized not in allowed:
            errors.append(f"{trid} : '{key}' invalide ({value!r}), valeurs possibles : {sorted(allowed)}")

    props = block.get("motif_properties")
    if not isinstance(props, dict):
        errors.append(f"{trid} : 'motif_properties' manquant")
        return
    for key in sorted(set(props) - MOTIF_PROPERTY_KEYS):
        warnings.append(f"{trid} : clé inconnue 'motif_properties.{key}' (ignorée)")

    groups = props.get("motif_groups")
    if not isinstance(groups, dict) or not groups:
        errors.append(f"{trid} : 'motif_properties.motif_groups' manquant ou vide")
        return
    for group_id, motifs in groups.items():
        if not isinstance(motifs, list) or not motifs:
            errors.append(f"{trid} [{group_id}] : groupe sans motif")
            continue
        for motif in motifs:
            if not isinstance(motif, str) or not _MOTIF_RE.match(motif):
                errors.append(f"{trid} [{group_id}] : motif invalide {motif!r}")

    thresholds = block.get("thresholds") or {}
    if not isinstance(thresholds, dict):
        errors.append(f"{trid} : 'thresholds' doit être un dictionnaire par groupe")
        thresholds = {}

    rules = block.get("structure_rules") or {}
    if isinstance(rules, list):
        rules = {}  # [] = aucune règle
    if not isinstance(rules, dict):
        errors.append(f"{trid} : 'structure_rules' doit être un dictionnaire par groupe ou []")
        rules = {}

    for group_id in sorted(set(thresholds) - set(groups)):
        errors.append(f"{trid} : thresholds définis pour un groupe inconnu '{group_id}'")
    for group_id in sorted(set(rules) - set(groups)):
        errors.append(f"{trid} : structure_rules définies pour un groupe inconnu '{group_id}'")

    structural = block.get("classification_mode") == "structural"
    if rules and not structural:
        warnings.append(f"{trid} : structure_rules ignorées (classification_mode = {block.get('classification_mode')!r})")

    for group_id in groups:
        where = f"{trid} [{group_id}]"
        group_thresholds = thresholds.get(group_id) or {}
        group_rules = rules.get(group_id) or []

        if not group_thresholds and not (structural and group_rules):
            errors.append(f"{where} : aucun seuil défini")
            continue

        labelled = {}
        for label, bounds in group_thresholds.items():
            if label not in label_priority:
                errors.append(f"{where} : label '{label}' absent de label_priority")
            checked = _check_range(bounds, f"{where} {label}", errors)
            if checked:
                labelled[label] = checked

        rule_ranges = []
        for i, rule in enumerate(group_rules, start=1):
            cond = rule.get("conditions") if isinstance(rule, dict) else None
            if not isinstance(cond, dict):
                errors.append(f"{where} règle {i} : 'conditions' manquant")
                continue
            label = cond.get("classification")
            if label not in label_priority:
                errors.append(f"{where} règle {i} : label '{label}' absent de label_priority")
            interruptions = cond.get("interruptions")
            if interruptions is not None and not isinstance(interruptions, bool):
                errors.append(f"{where} règle {i} : 'interruptions' doit être true, false ou absent")
            checked = _check_range(cond.get("repeat_range"), f"{where} règle {i}", errors)
            if checked:
                rule_ranges.append(checked)

        for l1, l2, lo, hi in overlaps(labelled):
            warnings.append(
                f"{where} : plages '{l1}' et '{l2}' qui se chevauchent sur {_format_gap(lo, hi)} "
                f"→ 'unclassified'"
            )

        covered = list(labelled.values()) + (rule_ranges if structural else [])
        for lo, hi in coverage_gaps(covered):
            warnings.append(f"{where} : valeurs {_format_gap(lo, hi)} hors de toute plage → 'unclassified'")


def validate_thresholds(data):
    """
    Valide le contenu chargé de clinical_thresholds.yaml.
    Retourne (errors, warnings), deux listes de messages.
    """
    errors, warnings = [], []

    if not isinstance(data, dict):
        return ["clinical_thresholds.yaml : contenu vide ou invalide"], warnings

    label_priority = _validate_label_priority(data, errors)

    low_depth = data.get("low_depth_threshold")
    if low_depth is not None and (not _is_int(low_depth) or low_depth < 0):
        errors.append(f"low_depth_threshold : entier positif ou vide attendu ({low_depth!r})")

    for trid, block in data.items():
        if trid in GLOBAL_KEYS:
            continue
        if not isinstance(block, dict):
            errors.append(f"{trid} : bloc de locus invalide")
            continue
        _validate_locus(trid, block, label_priority, errors, warnings)

    return errors, warnings
