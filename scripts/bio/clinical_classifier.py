from scripts.bio.labels import UNCLASSIFIED, NOTE_OUT_OF_RANGE, note_overlap


def in_range(value, bounds):
    """Bornes inclusives [min, max], None = non borné."""
    min_v, max_v = bounds
    if min_v is not None and value < min_v:
        return False
    if max_v is not None and value > max_v:
        return False
    return True


def clinical_group(data_group, clinical_group, repeat_mode, classification_mode):
    """
    Classification clinique d'un SEUL groupe.
    Retourne (label, note) ; note est None sauf pour 'unclassified'.

    Aucun repli silencieux : une valeur hors de toute plage, ou couverte par
    plusieurs plages de labels différents, donne 'unclassified'.
    """

    # 1) repeat_count clinique selon le repeat_mode du groupe
    if repeat_mode == "sum_with_interruptions":
        repeat_count = data_group.total_main_count_with
    else:
        repeat_count = data_group.total_main_count_without

    # 2) Mode structural → structure_rules d'abord (la première règle valide gagne)
    if classification_mode == "structural" and clinical_group.structure_rules:
        has_interruptions = (data_group.i_count > 0)

        for rule in clinical_group.structure_rules:
            cond = rule["conditions"]

            if not in_range(repeat_count, cond["repeat_range"]):
                continue

            # Vérification interruptions
            interruptions_required = cond.get("interruptions", None)
            if interruptions_required is not None and interruptions_required != has_interruptions:
                continue

            return cond["classification"], None

    # 3) Thresholds : toutes les plages sont examinées pour détecter les chevauchements
    matching = []
    for lbl, bounds in clinical_group.thresholds.items():
        if in_range(repeat_count, bounds) and lbl not in matching:
            matching.append(lbl)

    if len(matching) == 1:
        return matching[0], None
    if len(matching) > 1:
        return UNCLASSIFIED, note_overlap(matching)

    # 4) Aucune plage : non classable (jamais 'normal' par défaut)
    return UNCLASSIFIED, NOTE_OUT_OF_RANGE
