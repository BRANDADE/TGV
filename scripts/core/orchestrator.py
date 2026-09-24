import logging

from scripts.models.trid import TRID
from scripts.models.sample import Sample
from scripts.models.result import Result
from scripts.models.display import DisplayRow, DisplayDetails, DisplayExport, DisplayHtml


from scripts.core.result_builder import fill_raw_base, fill_clinical_base
from scripts.core.clinical_compute import genotype_value
from scripts.core.marking import mark_pathogenic_motifs, mark_pathogenic_segments, mark_pathogenic_repetition, mark_pathogenic_genotype
from scripts.core.rows import build_row_simple, build_row_clinical

from scripts.bio.motif_structure import decompose_repetition_without_interruptions, decompose_repetition_with_interruptions
from scripts.bio.clinical_classifier import clinical_group
from scripts.bio.labels import CALLED, NO_CALL, ABSENT, ABSENT_DISPLAY, UNCLASSIFIED, NOTE_MOTIF_DISCORDANCE, label_score

# Allèles appelés d'abord, puis non appelés, puis absents (locus haploïde)
_STATUS_RANK = {CALLED: 0, NO_CALL: 1, ABSENT: 2}


def to_int(value):
    """Entier ou None (valeur VCF manquante / vide)."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def depth_display(raw, status, threshold):
    """Profondeur affichée (avec ⚠ si sous le seuil) et profondeur brute exportée."""
    if status == ABSENT:
        return ABSENT_DISPLAY, ABSENT_DISPLAY
    value = to_int(raw)
    if value is None:
        return ".", "."
    shown = f"{value}"
    if threshold is not None and value < threshold:
        shown = f"\U000026A0 {value}"
    return shown, f"{value}"


def allele_sort_key(allele, clinical_cfg, min_label):
    """
    Clé de tri des allèles : (statut, génotype clinique affiché, AL).
    Locus sans configuration clinique (ou allèle sans génotype clinique) : AL seul.
    """
    inf = float("inf")
    clinical_value = inf
    if clinical_cfg is not None and allele.clinical is not None:
        clinical_value, _ = genotype_value(
            allele.clinical, clinical_cfg.genotype_display, clinical_cfg.pure_only, min_label
        )
    size = to_int(allele.size)
    return (_STATUS_RANK.get(allele.status, 3), clinical_value, size if size is not None else inf)


def order_alleles(sample, clinical_cfg, min_label):
    """
    Ordonne les allèles d'un sample : allèle 1 = plus petit génotype clinique affiché.
    Tri stable : à égalité, l'ordre TRGT est conservé. Les objets Allele entiers
    sont échangés (profondeur, pureté, méthylation et séquence suivent).
    """
    ordered = sorted(
        [sample.allele1, sample.allele2],
        key=lambda a: allele_sort_key(a, clinical_cfg, min_label),
    )
    sample.allele1, sample.allele2 = ordered


def raw_display(raw, status):
    """Valeur brute affichée pour un allèle ('-' si absent, '.' si manquante)."""
    if status == ABSENT:
        return ABSENT_DISPLAY
    return raw if raw not in (None, "") else "."


def process_result(analysis_input):
    label_priority = analysis_input.label_priority
    min_label = min(label_priority, key=label_priority.get)

    for trid_id, trid_global, sample in analysis_input.iter_items():

        if sample.result:
            continue

        # Ordre d'affichage des allèles fixé une fois pour toutes (UI, exports, IGV)
        order_alleles(sample, trid_global.clinical, min_label)
        a1 = sample.allele1
        a2 = sample.allele2

        # Un Result PAR TRID
        result = Result()

        # 1) RAW TRGT
        fill_raw_base(result, trid_id, trid_global, a1, a2)
        result.has_clinical = trid_global.clinical is not None

        # 2) Clinique si applicable
        if trid_global.clinical:
            fill_clinical_base(result, trid_global, a1, a2, min_label)

            # Motifs cliniques utilisés
            result.motifs_used1 = a1.clinical_motifs
            result.motifs_used2 = a2.clinical_motifs

        # 3) Lier le Result au Sample
        sample.result = result

    return None


def process_rows(analysis_input):
    """
    Construit les lignes simples TRGT pour l'UI à partir du DTO.
    """
    rows = []
    label_priority = analysis_input.label_priority
    min_label = min(label_priority, key=label_priority.get)

    # iter_items() retourne (trid_id, TRID_global, Sample)
    for trid_id, trid_global, sample in analysis_input.iter_items():
        
        a1 = sample.allele1
        a2 = sample.allele2

        clinical_cfg = trid_global.clinical
        if clinical_cfg is None:
            row = build_row_simple(
                trid_id=trid_id,
                trid=trid_global,
                a1=a1,
                a2=a2,
                paths=analysis_input.paths,
                sample_name=analysis_input.sample_name,
            )

        else:
           row = build_row_clinical(
                trid_id=trid_id,
                trid=trid_global,
                a1=a1,
                a2=a2,
                paths=analysis_input.paths,
                sample_name=analysis_input.sample_name,
                min_label=min_label
            )

        rows.append(row)

    return rows


def process_clinical(analysis_input):
    """
    Affiche la décomposition TRGT pour chaque TRID sélectionné
    et pour chaque allèle du sample.
    Remplit allele.trgt_groups[group_id] avec un TRGTGroupData structuré.

    Retourne la liste triée des TRID en discordance BED ↔ clinique
    (affichage à la charge de l'appelant : popup dans l'interface, log en CLI).
    """
    label_priority = analysis_input.label_priority

    logging.info(f"Executing clinical guidelines evaluation for sample: '{analysis_input.sample_name}'")

    # Liste des TRIDs en discordance BED ↔ clinique
    discordances = []

    for trid_id, trid_global, sample in analysis_input.iter_items():

        if sample.clinical_done:
            logging.debug(f"Clinical classification already computed for locus '{trid_id}'. Skipping recalculation.")
            continue

        clinical_cfg = trid_global.clinical
        if clinical_cfg is None:
            logging.debug(f"No clinical guidelines configuration found in YAML for locus '{trid_id}'. Skipping.")
            continue

        repeat_mode = clinical_cfg.repeat_mode
        classification_mode = clinical_cfg.classification_mode
        motif_groups = clinical_cfg.groups

        alleles = [sample.allele1, sample.allele2]
        sample.clinical_done = True

        for idx, allele in enumerate(alleles, start=1):

            # Allèle non appelé (no-call) ou absent (locus haploïde) : rien à classer
            if not allele.is_called:
                logging.debug(f"Allele {idx} of locus '{trid_id}' is '{allele.status}'. No classification.")
                continue

            logging.debug(f"Processing allele {idx} for locus '{trid_id}':")
            logging.debug(f"  Repetitions: {allele.sequence.repetitions}")
            logging.debug(f"  Segmentation: {allele.sequence.segmentation}")
            logging.debug(f"  Sequence preview: {allele.sequence.sequence[:50]}...")

            allele.trgt_groups = {}
            allele.clinical = None
            allele.clinical_label = None
            allele.clinical_note = None

            best_group = None
            best_score = -1

            for group_id, group in motif_groups.items():

                # Décomposition TRGT
                if repeat_mode == "sum_with_interruptions":
                    data = decompose_repetition_with_interruptions(
                        repetitions=allele.sequence.repetitions,
                        segmentation=allele.sequence.segmentation_complete,
                        groups=group.motifs
                    )
                else:
                    data = decompose_repetition_without_interruptions(
                        repetitions=allele.sequence.repetitions,
                        groups=group.motifs
                    )

                # --- CAS 1 : data est None → mismatch motifs BED / clinique ---
                if data is None:
                    logging.warning(f"Locus '{trid_id}' mismatch: No motifs of clinical group '{group_id}' found in TRGT outputs.")
                    allele.trgt_groups[group_id] = None
                    continue

                # Classification clinique (jamais de repli silencieux sur 'normal')
                clinical_label, clinical_note = clinical_group(
                    data_group=data,
                    clinical_group=group,
                    repeat_mode=repeat_mode,
                    classification_mode=classification_mode,
                )

                data.clinical = clinical_label
                data.clinical_note = clinical_note

                # Debug
                logging.debug(f"  Clinical group evaluation details - Group: {group_id} | Motifs: {group.motifs}")
                logging.debug(f"  Decomposition: {data}")

                allele.trgt_groups[group_id] = data

                # Score clinique
                score = label_score(clinical_label, label_priority)
                
                # Total repeats clinique
                if data.total_main_count_with is not None:
                    total_repeats = data.total_main_count_with
                else:
                    total_repeats = data.total_main_count_without

                # Sélection du meilleur groupe
                if best_group is None:
                    best_group = group_id
                    best_score = score
                    best_repeats = total_repeats
                else:
                    if score > best_score:
                        best_group = group_id
                        best_score = score
                        best_repeats = total_repeats
                    elif score == best_score and total_repeats > best_repeats:
                        best_group = group_id
                        best_repeats = total_repeats

            # --- Aucun groupe valide trouvé : allèle conservé, marqué non classable ---
            if best_group is None:
                logging.error(f"No valid clinical threshold group resolved for allele {idx} of locus '{trid_id}'.")
                allele.clinical_label = UNCLASSIFIED
                allele.clinical_note = NOTE_MOTIF_DISCORDANCE
                discordances.append(trid_id)
                continue

            # Stockage du groupe gagnant
            allele.clinical = allele.trgt_groups[best_group]
            allele.clinical_label = allele.clinical.clinical
            allele.clinical_note = allele.clinical.clinical_note
            allele.clinical_motifs = motif_groups[best_group].motifs
            if allele.clinical_label == UNCLASSIFIED:
                logging.warning(f"Allele {idx} of locus '{trid_id}' is unclassified: {allele.clinical_note}.")

    discordances = sorted(set(discordances))
    if discordances:
        logging.warning(f"Clinical/genomic discordances resolved for loci: {discordances}")

    return discordances


def process_display(result, clinical_cfg, low_depth_threshold):
    # Si déjà calculé → ne rien faire
    if result.display_row is not None:
        return

    valid_threshold = None  # Default is None (no warning marker displayed)

    if low_depth_threshold is not None:
        try:
            valid_threshold = int(low_depth_threshold)
        except (ValueError, TypeError) as e:
            logging.warning(
                f"The configuration parameter 'low_depth_threshold' ('{low_depth_threshold}') "
                f"is not a valid integer ({e}). Low coverage visual warning markers are disabled."
            )
            valid_threshold = None

    pathogenic = None
    if clinical_cfg:
        pathogenic = clinical_cfg.pathogenic_motifs

    result.display_row = DisplayRow()
    result.display_details = DisplayDetails()
    result.display_export = DisplayExport()
    result.display_html = DisplayHtml()

    # --- Locus ---
    result.display_row.locus = result.locus
    result.display_details.locus = result.locus
    result.display_export.locus = result.locus
    result.display_html.locus = result.locus

    # --- Profondeur ---
    depth1, result.display_export.depth1 = depth_display(result.depth1_raw, result.status1, valid_threshold)
    depth2, result.display_export.depth2 = depth_display(result.depth2_raw, result.status2, valid_threshold)
    result.display_row.depth = f"{depth1} / {depth2}"
    result.display_details.depth = f"{depth1} / {depth2}"
    result.display_html.depth = f"{depth1} / {depth2}"

    # --- Taille ---
    size1 = raw_display(result.size1_raw, result.status1)
    size2 = raw_display(result.size2_raw, result.status2)
    size = f"{size1} / {size2}"
    result.display_row.size = size
    size_extended = (
        f"{size1} ({raw_display(result.range_size1_raw, result.status1)}) / "
        f"{size2} ({raw_display(result.range_size2_raw, result.status2)})"
    )
    result.display_details.size = size_extended

    # --- Motifs TRGT ---
    # 1) Export CIL → brut, propre
    motifs_str = ", ".join(result.motifs_raw)
    result.display_export.motifs = motifs_str

    # 2) UI / HTML → marquage si clinique
    if pathogenic:
        motifs_str_ui = mark_pathogenic_motifs(result.motifs_raw, pathogenic, ui=True)
        motifs_str_html = mark_pathogenic_motifs(result.motifs_raw, pathogenic)
    else:
        motifs_str_ui = motifs_str
        motifs_str_html = motifs_str
        
    # 3) Affectation UI
    result.display_row.motifs = motifs_str_ui
    result.display_details.motifs = motifs_str_ui

    # 4) Affectation HTML
    result.display_html.motifs = motifs_str_html


    # --- Motifs utilisés (cliniques) ---
    motifs_use1 = result.motifs_used1 or []
    motifs_use2 = result.motifs_used2 or []

    if clinical_cfg:
        motifs_use1_ui = mark_pathogenic_motifs(motifs_use1, pathogenic, ui=True)
        motifs_use2_ui = mark_pathogenic_motifs(motifs_use2, pathogenic, ui=True)
    else:
        motifs_use1_ui = ", ".join(motifs_use1)
        motifs_use2_ui = ", ".join(motifs_use2)

    result.display_details.motifs_use1 = motifs_use1_ui
    result.display_details.motifs_use2 = motifs_use2_ui

    # --- Genotype ----
    if result.genotype1_raw or result.genotype2_raw:
        genotype = f"{result.genotype1_raw} / {result.genotype2_raw}"
    else:
        genotype = ""

    # Export = brut
    result.display_export.genotype = genotype

    # UI / HTML marqués si clinique
    if pathogenic:
        genotype_ui = mark_pathogenic_genotype(genotype, pathogenic, ui=True)
        genotype_html = mark_pathogenic_genotype(genotype, pathogenic, ui=False)
    else:
        genotype_ui = genotype
        genotype_html = genotype

    # Affectation UI
    result.display_row.genotype = genotype_ui
    result.display_details.genotype = genotype_ui

    # Affectation HTML
    result.display_html.genotype = genotype_html

    # --- Classification ---

    # Allele 1
    if result.classification1_raw:
        classif1 = result.classification1_bio if result.classification1_bio else result.classification1_raw
    else:
        classif1 = ""

    # Allele 2
    if result.classification2_raw:
        classif2 = result.classification2_bio if result.classification2_bio else result.classification2_raw
    else:
        classif2 = ""

    # Format final
    classification = f"{classif1} / {classif2}" if (classif1 or classif2) else ""

    # Affectation UI
    result.display_row.classification = classification
    result.display_details.classification = classification

    # Export = brut (pas de marquage)
    result.display_export.classification = classification

    # HTML = brut ou marqué si tu veux plus tard
    result.display_html.classification = classification

    # --- Répétitions ---

    # Allele 1 : choisir la bonne source
    if result.rep1_clinical_raw:
        rep1_base = result.rep1_clinical_raw
    else:
        rep1_base = result.rep1_raw or ""

    # Allele 2 : choisir la bonne source
    if result.rep2_clinical_raw:
        rep2_base = result.rep2_clinical_raw
    else:
        rep2_base = result.rep2_raw or ""

    # Export = brut, jamais marqué
    result.display_export.rep1 = rep1_base
    result.display_export.rep2 = rep2_base

    # UI / HTML marqués si clinique
    if pathogenic:
        rep1_ui = mark_pathogenic_repetition(rep1_base, pathogenic, ui=True)
        rep2_ui = mark_pathogenic_repetition(rep2_base, pathogenic, ui=True)

        rep1_html = mark_pathogenic_repetition(rep1_base, pathogenic, ui=False)
        rep2_html = mark_pathogenic_repetition(rep2_base, pathogenic, ui=False)
    else:
        rep1_ui = rep1_base
        rep2_ui = rep2_base
        rep1_html = rep1_base
        rep2_html = rep2_base

    # Affectation UI
    result.display_row.rep1 = rep1_ui
    result.display_row.rep2 = rep2_ui
    result.display_details.rep1 = rep1_ui
    result.display_details.rep2 = rep2_ui

    # Affectation HTML
    result.display_html.rep1 = rep1_html
    result.display_html.rep2 = rep2_html

    # --- Segmentation ---
    seg1_base = result.seg1_raw or ""
    seg2_base = result.seg2_raw or ""

    # Export = brut
    result.display_export.seg1 = seg1_base
    result.display_export.seg2 = seg2_base

    # UI / HTML marqués si clinique
    if clinical_cfg:
        seg1_ui = mark_pathogenic_segments(seg1_base, pathogenic, ui=True)
        seg2_ui = mark_pathogenic_segments(seg2_base, pathogenic, ui=True)

        seg1_html = mark_pathogenic_segments(seg1_base, pathogenic, ui=False)
        seg2_html = mark_pathogenic_segments(seg2_base, pathogenic, ui=False)
    else:
        seg1_ui = seg1_base
        seg2_ui = seg2_base
        seg1_html = seg1_base
        seg2_html = seg2_base

    # Affectation UI
    result.display_row.seg1 = seg1_ui
    result.display_row.seg2 = seg2_ui
    result.display_details.seg1 = seg1_ui
    result.display_details.seg2 = seg2_ui

    # Affectation HTML
    result.display_html.seg1 = seg1_html
    result.display_html.seg2 = seg2_html

    # --- Interruptions ---- # 
    result.display_details.interruptions1 = result.inter1_raw or ""
    result.display_details.interruptions2 = result.inter2_raw or ""


    # --- Sequence ---- # 
    result.display_details.sequence1 = result.sequence1_raw or ""
    result.display_details.sequence2 = result.sequence2_raw or ""


    # --- Purity --- #
    purity = f"{result.purity1} / {result.purity2}"
    result.display_details.purity = purity

    # --- Methylation --- #
    methylation = f"{result.methylation1} / {result.methylation2}"
    result.display_details.methylation = methylation