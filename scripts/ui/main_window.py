import PySimpleGUI as sg
import os
import sys
import yaml
import re
import logging

from scripts.models.run import Run

from scripts.core.vcf_loader import list_vcfs
from scripts.core.vcf_parser import parse_vcf_for_sample
from scripts.core.utils import get_analysis_prefix
from scripts.core.artifact_lookup import sample_id_from_vcf_name
from scripts.core.config_manager import get_safe_config_path

from scripts.bio.clinical_thresholds_loader import load_clinical_thresholds, load_trid_aliases
from scripts.bio.clinical_config_validator import validate_thresholds
from scripts.core.analysis import build_run_trids, resolve_panel, run_analysis
from scripts.core.provenance import read_vcf_header, thresholds_sha256, tgv_provenance

from scripts.ui.results_window import show_results_window
from scripts.ui.igv import is_online, get_asset_path
from scripts.ui.make_report import open_report_on_the_fly

# ---------------------------------------------------------
# Chargement robuste du fichier buttons_panel.yaml
# ---------------------------------------------------------
def load_button_panels():
    # Appel de la fonction de chemin sécurisée
    path = get_safe_config_path("buttons_panel.yaml")

    if not os.path.exists(path):
        logging.info(f"No 'buttons_panel.yaml' found at {path}. No panel loaded.")
        return {}

    try:
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        return data if isinstance(data, dict) else {}
    except Exception as e:
        logging.warning(f"Failed to load buttons_panel.yaml: {e}")
        return {}


# ---------------------------------------------------------
# Flow layout automatique
# ---------------------------------------------------------
def build_panel_rows(panel_names, window_width=600, button_padding=40):
    rows, row, current_width = [], [], 0

    for name in panel_names:
        estimated_width = len(name) * 8 + button_padding
        if current_width + estimated_width > window_width:
            rows.append(row)
            row, current_width = [], 0
        row.append(sg.Button(name, key=f"-PANEL-{name}-"))
        current_width += estimated_width

    if row:
        rows.append(row)

    return rows


# ---------------------------------------------------------
# Mise à jour TRIDs sélectionnés
# ---------------------------------------------------------
def update_trid_selected(window):
    selected = window.metadata.get("selected_trids", [])
    readable_map = window.metadata.get("readable_to_trid", {})
    trid_to_readable = {v: k for k, v in readable_map.items()}
    readable_list = [trid_to_readable[t] for t in selected]
    window["-TRID-SELECTED-"].update(values=readable_list)


# ---------------------------------------------------------
# Récupère les TRIDs ayant un bloc clinique
# ---------------------------------------------------------
def get_trids_with_clinical(window):
    run = window.metadata.get("run")
    if not run:
        return set()
    return {
        trid_id
        for trid_id, trid in run.trids.items()
        if getattr(trid, "clinical", None) is not None
    }


# ---------------------------------------------------------
# FENÊTRE PRINCIPALE
# ---------------------------------------------------------
def run_main_window():
    logging.info("Opening main window interface")

    sg.theme("SystemDefault")
    # Ligne de titre stylisée à insérer au début de votre layout
    title_row = [
        # T - RGT
        sg.Push(),
        sg.Text("T", font=("Helvetica", 24, "bold"), text_color="#c8336a", pad=(0, 0)),
        sg.Text("RGT", font=("Helvetica", 14, "bold"), text_color="#5a4050", pad=((0, 15), 0)),
        
        # G - lobal
        sg.Text("G", font=("Helvetica", 24, "bold"), text_color="#c8336a", pad=(0, 0)),
        sg.Text("lobal", font=("Helvetica", 14, "bold"), text_color="#5a4050", pad=((0, 15), 0)),
        
        # V - iewer
        sg.Text("V", font=("Helvetica", 24, "bold"), text_color="#c8336a", pad=(0, 0)),
        sg.Text("iewer", font=("Helvetica", 14, "bold"), text_color="#5a4050", pad=(0, 0)),
        sg.Push(),
]

    layout = [
        title_row,
        [sg.Text("Sélection du fichier TRGT (trgt_vcfs.zip)")],
        [sg.Input(key="-ZIP-", enable_events=True), sg.FileBrowse("Parcourir")],
        [sg.Button("Rapport QC", key="-QC-REPORT-", disabled=True)],

        [sg.Text("Sélection du génome de référence (.fa / .fasta) [Optionnel]")],
        [sg.Input(key="-FASTA-", enable_events=True), sg.FileBrowse("Parcourir")],
        
        [sg.Text("Patient à analyser")],
        [sg.Input(key="-SEARCH-", enable_events=True, size=(40, 1))],
        [sg.Combo([], key="-SAMPLE-", size=(40, 1), readonly=True, enable_events=True)],

        [sg.Text("Recherche locus")],
        [sg.Input(key="-SEARCH-TRID-", enable_events=True, size=(40, 1))],
        [
            sg.Combo([], key="-TRID-COMBO-", size=(40, 1), readonly=True),
            sg.Button("Ajouter", key="-ADD-TRID-"),
        ],

        [sg.Frame(
            "Locus sélectionnés",
            [
                [sg.Listbox([], key="-TRID-SELECTED-", size=(40, 6),
                            select_mode=sg.SELECT_MODE_EXTENDED)],
                [sg.Button("Retirer", key="-DEL-TRID-", size=(12, 1))],
            ],
        )],

        [sg.Frame(
            "Actions sur les locus",
            [
                [sg.Button("Tout cocher", key="-SELECT-ALL-", size=(15, 1)),
                 sg.Button("Tout décocher", key="-UNSELECT-ALL-", size=(15, 1))]
            ],
        )],
    ]

    # Panels YAML
    panels = load_button_panels()
    if panels:
        panel_rows = build_panel_rows(list(panels.keys()))
        layout.append([sg.Frame("Panels", panel_rows)])

    logo_path = get_asset_path("logo_chu_Nimes.png")

    layout.extend([
        [
            sg.Column([
                [sg.Button("Lancer l'analyse")],
                [sg.Text("", key="-STATUS-", text_color="red", font=("Helvetica", 9), pad=(0,0))]
            ], vertical_alignment="bottom", pad=(0,0)), 

            sg.Push(),
            sg.Column([
                
                [sg.Image(filename=logo_path, subsample=3) if os.path.exists(logo_path) else sg.Text("")]
            ], vertical_alignment="bottom", element_justification="right", pad=(0,0)) 
        ],
        [
            sg.Text("by Corentin Marco", font=("Helvetica", 8), text_color="gray", pad=((5,0),(0,0)))
        ],
    ])

    # Modification du titre officiel de la fenêtre
    window = sg.Window("TGV - TRGT Global Viewer", layout, finalize=True)
    

    window.metadata = {
        "all_samples": [],
        "all_trids": [],
        "readable_to_trid": {},
        "selected_trids": [],
        "button_panels": panels
    }

    # ---------------------------------------------------------
    # BOUCLE D'ÉVÉNEMENTS
    # ---------------------------------------------------------
    # Étape 1 : Tester la connexion une seule fois au démarrage de la fenêtre
    logging.info("Checking Internet connection for IGV dependency...")
    online_status = is_online()
    logging.info(f"Internet connectivity status: {'Online' if online_status else 'Offline'}")

    fasta_path = None
    while True:
        event, values = window.read()

        if event == sg.WINDOW_CLOSED:
            logging.info("Closing main window interface")

            break

        if event == "-FASTA-":
            fasta_path = values["-FASTA-"]
            if not fasta_path:
                continue

            # Vérifier extension
            logging.info(f"Reference genome FASTA file selected: {fasta_path}")
            fasta_path_lower = fasta_path.lower()
            if not (fasta_path_lower.endswith(".fasta") or fasta_path_lower.endswith(".fa")):
                logging.warning(f"Invalid FASTA file format: {fasta_path}")
                sg.popup_error("Le fichier doit être au format .fasta ou .fa")
                continue

            # Vérifier existence du fichier FASTA
            if not os.path.isfile(fasta_path):
                logging.warning(f"FASTA file not found on disk: {fasta_path}")
                sg.popup_error("Le fichier FASTA n'existe pas.")
                continue

            # Vérifier existence du fichier .fai
            fai_path = fasta_path + ".fai"
            if not os.path.isfile(fai_path):
                logging.warning(f"Missing FASTA index (.fai) for: {fasta_path}")
                sg.popup_error(
                    f"Le fichier d'index FASTA (.fai) est manquant :\n{fai_path}\n\n"
                    "Génère-le avec :\n\nsamtools faidx fichier.fasta"
                )
                continue

        # ---------------------------------------------------------
        # Chargement ZIP
        # ---------------------------------------------------------
        if event == "-ZIP-":
            zip_path = values["-ZIP-"]
            if not zip_path:
                continue

            logging.info(f"Loading main TRGT ZIP archive: {zip_path}")

            # Nouveau run : rien de l'ancien run ne doit subsister (sélection de loci,
            # patient, allèles parsés)
            window.metadata.pop("run", None)
            window.metadata["current_trids"] = {}
            window.metadata["current_sample"] = None
            window.metadata["selected_trids"] = []
            window.metadata["all_samples"] = []
            window["-SAMPLE-"].update(values=[], value="")
            update_trid_selected(window)
            window["-STATUS-"].update("")

            vcfs = list_vcfs(zip_path)
            if not vcfs:
                logging.error(f"Selected ZIP archive contains no .trgt.vcf files: {zip_path}")
                sg.popup_error("Ce ZIP ne contient aucun fichier .trgt.vcf.")
                continue

            run_name = get_analysis_prefix(zip_path)
            logging.info(f"Analysis run name prefix identified: {run_name}")
            run = Run(name=run_name, vcf_zip=zip_path)

            # ---------------------------------------------------------
            # Détection ZIP associés
            # ---------------------------------------------------------

            zip_filename = os.path.basename(zip_path)

            # On retire uniquement le suffixe final
            # → "analysis-S18-30.04.2026-"
            prefix = re.sub(r"(?i)trgt_vcfs\.zip", "", zip_filename)

            # Répertoire contenant les fichiers
            base_dir = os.path.dirname(zip_path)

            # Construction explicite des chemins
            run.spanning_bam_zip      = os.path.join(base_dir, prefix + "spanning_BAM.zip")
            run.repeat_reads_zip      = os.path.join(base_dir, prefix + "repeat_reads.zip")
            run.motifs_allele_zip     = os.path.join(base_dir, prefix + "trgt_motifs_allele.zip")
            run.motifs_waterfall_zip  = os.path.join(base_dir, prefix + "trgt_motifs_waterfall.zip")
            run.meth_allele_zip       = os.path.join(base_dir, prefix + "trgt_meth_allele.zip")
            run.meth_waterfall_zip    = os.path.join(base_dir, prefix + "trgt_meth_waterfall.zip")
            run.qc_zip                = os.path.join(base_dir, prefix + "QC.zip")

            # Vérification d'existence
            if not os.path.exists(run.spanning_bam_zip):
                run.spanning_bam_zip = None
                logging.info("Associated Spanning BAM archive: NOT FOUND")
            else :
                logging.info(f"Associated Spanning BAM archive found: {run.spanning_bam_zip}")

            # Vérification d'existence
            if not os.path.exists(run.repeat_reads_zip):
                run.repeat_reads_zip = None
                logging.info("Associated Repeat Reads archive: NOT FOUND")
            else:
                logging.info(f"Associated Repeat Reads archive found: {run.repeat_reads_zip}")

            if not os.path.exists(run.motifs_allele_zip):
                run.motifs_allele_zip = None
                logging.info("Associated Motifs Allele archive: NOT FOUND")
            else:
                logging.info(f"Associated Motifs Allele archive found: {run.motifs_allele_zip}")

            if not os.path.exists(run.motifs_waterfall_zip):
                run.motifs_waterfall_zip = None
                logging.info("Associated Motifs Waterfall archive: NOT FOUND")
            else:
                logging.info(f"Associated Motifs Waterfall archive found: {run.motifs_waterfall_zip}")

            if not os.path.exists(run.meth_allele_zip):
                run.meth_allele_zip = None
                logging.info("Associated Methylation Allele archive: NOT FOUND")
            else:
                logging.info(f"Associated Methylation Allele archive found: {run.meth_allele_zip}")

            if not os.path.exists(run.meth_waterfall_zip):
                run.meth_waterfall_zip = None
                logging.info("Associated Methylation Waterfall archive: NOT FOUND")
            else:
                logging.info(f"Associated Methylation Waterfall archive found: {run.meth_waterfall_zip}")

            if run.qc_zip and os.path.exists(run.qc_zip):
                logging.info(f"Associated Quality Control (QC) archive found: {run.qc_zip}")
                window["-QC-REPORT-"].update(disabled=False)
            else:
                logging.info("Associated Quality Control (QC) archive: NOT FOUND")
                run.qc_zip = None
                window["-QC-REPORT-"].update(disabled=True)

            # Construire la liste des noms de fichiers VCF
            display = [os.path.basename(v) for v in vcfs]

            # Mettre à jour la combo avec une LISTE
            window["-SAMPLE-"].update(values=display)

            # Stocker la liste complète
            window.metadata["all_samples"] = display

            # ---------------------------------------------------------
            # Configuration clinique : validation avant tout usage
            # ---------------------------------------------------------
            thresholds_data = load_clinical_thresholds()
            errors, warnings = validate_thresholds(thresholds_data)

            for w in warnings:
                logging.warning(f"clinical_thresholds.yaml: {w}")

            if errors:
                for e in errors:
                    logging.error(f"clinical_thresholds.yaml: {e}")
                sg.popup_error(
                    "clinical_thresholds.yaml est invalide ; le run n'est pas chargé.\n\n"
                    + "\n".join(f" - {e}" for e in errors),
                    title="Configuration clinique invalide",
                    keep_on_top=True,
                )
                window["-SAMPLE-"].update(values=[], value="")
                window.metadata["all_samples"] = []
                continue

            if warnings and not window.metadata.get("yaml_warnings_shown"):
                window.metadata["yaml_warnings_shown"] = True
                sg.popup_scrolled(
                    "Avertissements de configuration (clinical_thresholds.yaml).\n"
                    "Les valeurs concernées seront classées 'unclassified'.\n\n"
                    + "\n".join(f" - {w}" for w in warnings),
                    title="Configuration clinique",
                    size=(100, 20),
                )

            label_priority = thresholds_data["label_priority"]
            low_depth_threshold = thresholds_data.get("low_depth_threshold", None)
            logging.info(f"clinical_thresholds.yaml SHA-256: {thresholds_sha256()}")

            # ---------------------------------------------------------
            # Création des TRIDs globaux + remplissage infos immuables + clinique
            # ---------------------------------------------------------
            trids, diseases, run.trids = build_run_trids(zip_path, thresholds_data, load_trid_aliases())
            logging.info(f"Initialized {len(trids)} genomic loci structures.")

            for trid_id, t in run.trids.items():
                logging.debug(f"Locus: {trid_id} | Chrom: {t.chrom} | Coords: {t.start}-{t.end} | Motifs: {t.motifs}")
                if t.clinical:
                    logging.debug(f"  Clinical Config: Mode={t.clinical.classification_mode} | Groups={list(t.clinical.groups.keys())}")

            logging.info("Successfully initialized all genomic loci structures.")

            # ---------------------------------------------------------
            # Stockage global
            # ---------------------------------------------------------
            window.metadata["run"] = run
            window.metadata["all_trids"] = trids
            window.metadata["readable_to_trid"] = diseases
            window.metadata["label_priority"] = label_priority

            window["-TRID-COMBO-"].update(values=sorted(diseases.keys()))

        # --------------------------------------------------------- 
        # Événement : Ouvrir le Rapport QC
        # ---------------------------------------------------------
        if event == "-QC-REPORT-":
            run = window.metadata.get("run")
            if run and run.qc_zip:
                open_report_on_the_fly(run.qc_zip)

        # ---------------------------------------------------------
        # Recherche patient
        # ---------------------------------------------------------
        
        if event == "-SEARCH-":
            query = values["-SEARCH-"].lower()
            all_samples = window.metadata.get("all_samples", [])

            # Filtrage
            filtered = [s for s in all_samples if query in s.lower()] if query else all_samples
            window["-SAMPLE-"].update(values=filtered)

            # Auto‑sélection si un seul résultat
            if len(filtered) == 1:
                sample = filtered[0]

                # Sélectionner SANS écraser la liste
                window["-SAMPLE-"].update(value=sample)

                # Déclencher le parsing
                window.write_event_value("-SAMPLE-", sample)

                # Effacer la recherche
                window["-SEARCH-"].update("")

                # Rétablir la liste complète
                window["-SAMPLE-"].update(values=all_samples, value=sample)

               #Retire le focus 
                window["-SAMPLE-"].set_focus()


        # ---------------------------------------------------------
        # Sélection patient
        # ---------------------------------------------------------
        if event == "-SAMPLE-":
            sample_name = values.get("-SAMPLE-")
            if not sample_name:
                logging.warning("User triggered patient selection but no sample name was found")
                continue

            run = window.metadata.get("run")
            if not run:
                continue

            # Parser le VCF → retourne { trid_id : Sample }
            sample_trids = parse_vcf_for_sample(
                zip_path=run.vcf_zip,
                vcf_filename=sample_name,
                global_trids=run.trids,
            )

            # Attacher chaque Sample au TRID correspondant
            for trid_id, sample_obj in sample_trids.items():
                if trid_id in run.trids:
                    run.trids[trid_id].samples[sample_name] = sample_obj

            # Stockage UI
            window.metadata["current_trids"] = sample_trids
            window.metadata["current_sample"] = sample_name

            # Mise à jour TRID lisibles
            readable_map = window.metadata["readable_to_trid"]
            readable_list = sorted(readable_map.keys())
            window["-TRID-COMBO-"].update(values=readable_list)

            # Log clinique de traçabilité
            logging.info(f"Patient sample selected: '{sample_name}' - Extracted {len(sample_trids)} active loci.")

            # --- DEBUG ---
            logging.debug(f"Detailed sequences for sample: {sample_name}")
            for trid_id, sample_obj in sample_trids.items():
                a1, a2 = sample_obj.allele1, sample_obj.allele2
                logging.debug(f"  Locus: {trid_id}")
                logging.debug(f"    Allele 1: Size={a1.size} | Consensus={a1.sequence.sequence} | MC={a1.sequence.repetitions}")
                logging.debug(f"    Allele 2: Size={a2.size} | Consensus={a2.sequence.sequence} | MC={a2.sequence.repetitions}")


        # ---------------------------------------------------------
        # Recherche TRID lisible
        # ---------------------------------------------------------
        if event == "-SEARCH-TRID-":
            query = values["-SEARCH-TRID-"].lower()
            readable_map = window.metadata.get("readable_to_trid", {})
            readable_list = sorted(readable_map.keys())

            filtered = [r for r in readable_list if query in r.lower()] if query else readable_list
            window["-TRID-COMBO-"].update(values=filtered)

            # Auto‑sélection si un seul TRID trouvé
            if len(filtered) == 1:
                readable = filtered[0]

                # Sélectionner SANS écraser la liste
                window["-TRID-COMBO-"].update(values=readable_list,value=readable)

                # Effacer la recherche
                window["-SEARCH-TRID-"].update("")

               #  Retire le focus
                window["-TRID-COMBO-"].set_focus()

        # ---------------------------------------------------------
        # Ajouter TRID
        # ---------------------------------------------------------
        if event == "-ADD-TRID-":
            readable = values["-TRID-COMBO-"]
            if readable:
                trid = window.metadata["readable_to_trid"][readable]
                if trid not in window.metadata["selected_trids"]:
                    window.metadata["selected_trids"].append(trid)
                    logging.info(f"Locus added to active list: {readable} ({trid})")
                update_trid_selected(window)

            # Effacer la sélection dans la combo TRID
            window["-TRID-COMBO-"].update(value="")

            # Effacer aussi la recherche TRID
            window["-SEARCH-TRID-"].update("")


        # ---------------------------------------------------------
        # Tout cocher / décocher
        # ---------------------------------------------------------
        if event == "-SELECT-ALL-":
            readable_map = window.metadata["readable_to_trid"]
            all_trids = list(readable_map.values())

            # 1. TRIDs avec clinique
            trids_clinical = get_trids_with_clinical(window)
            first = [t for t in all_trids if t in trids_clinical]

            # 2. TRIDs sans clinique
            second = [t for t in all_trids if t not in trids_clinical]

            # Ordre final
            window.metadata["selected_trids"] = first + second

            update_trid_selected(window)


        if event == "-UNSELECT-ALL-":
            window.metadata["selected_trids"] = []
            update_trid_selected(window)

        # ---------------------------------------------------------
        # Retirer TRID
        # ---------------------------------------------------------
        if event == "-DEL-TRID-":
            readable_selected = values["-TRID-SELECTED-"] or []
            readable_map = window.metadata["readable_to_trid"]
            selected = window.metadata["selected_trids"]
            for readable in readable_selected:
                trid = readable_map[readable]
                if trid in selected:
                    selected.remove(trid)
                    logging.info(f"Locus removed from active list: {readable} ({trid})")
            update_trid_selected(window)

        # ---------------------------------------------------------
        # Panels YAML
        # ---------------------------------------------------------
        if event.startswith("-PANEL-"):
            panel_name = event.replace("-PANEL-", "").replace("-", "")
            panels = window.metadata["button_panels"]
            if panel_name in panels:
                requested = panels[panel_name]
                run = window.metadata.get("run")
                existing, missing = resolve_panel(requested, run.trids if run else {})

                window.metadata["selected_trids"] = existing
                update_trid_selected(window)

                if missing:
                    readable_map = window.metadata["readable_to_trid"]
                    missing_readable = []
                    for trid in missing:
                        readable = next((k for k, v in readable_map.items() if v == trid), trid)
                        missing_readable.append(readable)
                    sg.popup("Locus absents du VCF :\n" + "\n".join(missing_readable), keep_on_top=True)

        # ---------------------------------------------------------
        # Lancer l'analyse
        # ---------------------------------------------------------
        if event == "Lancer l'analyse":
            run = window.metadata.get("run")
            sample_name = values["-SAMPLE-"]
            selected_trids = window.metadata.get("selected_trids")

            if not run or not sample_name:
                logging.warning("User tried to launch analysis but no ZIP archive or patient sample was selected.")
                window["-STATUS-"].update("Veuillez sélectionner un ZIP et un patient", text_color="red")
                continue

            if not selected_trids:
                logging.warning("User tried to launch analysis but no loci were selected.")
                window["-STATUS-"].update("Veuillez sélectionner au moins un locus", text_color="red")
                continue

            # Seuls les loci présents chez ce patient sont analysés ; les autres sont signalés.
            # Les allèles doivent être ceux du patient affiché (re-parsing si besoin).
            if window.metadata.get("current_sample") != sample_name:
                window.metadata["current_trids"] = parse_vcf_for_sample(run.vcf_zip, sample_name, run.trids)
                window.metadata["current_sample"] = sample_name
            sample_trids = window.metadata["current_trids"]
            paths = {
                "vcf": run.vcf_zip,
                "repeat_reads": run.repeat_reads_zip,
                "spanning_bam": run.spanning_bam_zip,
                "motifs_allele": run.motifs_allele_zip,
                "motifs_waterfall": run.motifs_waterfall_zip,
                "meth_allele": run.meth_allele_zip,
                "meth_waterfall": run.meth_waterfall_zip,
                "genome_fasta": fasta_path,
            }

            logging.info(f"Initiating clinical analysis for patient: '{sample_name}' - Loci selected: {len(selected_trids)}")
            for path_key, path_val in paths.items():
                logging.debug(f"    Param path -> {path_key}: {path_val}")

            label_priority = window.metadata["label_priority"]
            results, discordances, missing = run_analysis(
                sample_name, sample_trids, run.trids, selected_trids,
                label_priority, low_depth_threshold, paths=paths,
            )

            if discordances:
                sg.popup_error(
                    "Discordance entre les motifs TRGT (BED) et les seuils cliniques définis dans "
                    "clinical_thresholds.yaml pour les locus suivants :\n\n"
                    + "\n".join(f" - {trid}" for trid in discordances),
                    title="Discordance clinique détectée",
                )

            if missing:
                trid_to_readable = {v: k for k, v in window.metadata["readable_to_trid"].items()}
                missing_readable = [trid_to_readable.get(t, t) for t in missing]
                logging.warning(f"Loci absent from the VCF of '{sample_name}': {missing}")
                window["-STATUS-"].update(
                    "Absents du VCF du patient : " + ", ".join(missing_readable), text_color="red"
                )
            else:
                window["-STATUS-"].update("")

            if not results:
                continue

            # Appelle l'UI avec la liste de Result
            logging.info(f"Analysis completed. Successfully processed {len(results)}/{len(selected_trids)} loci.")
            logging.info(f"Opening results visualization window for patient: '{sample_name}'")

            run_id = run.name if run else ""

            # Provenance affichée dans l'export : TGV, TRGT (en-tête VCF), catalogue, YAML
            provenance = tgv_provenance()
            try:
                provenance.update(read_vcf_header(run.vcf_zip, sample_name))
            except Exception as e:
                logging.warning(f"Cannot read TRGT header of '{sample_name}': {e}")

            show_results_window(
                sample_name=sample_id_from_vcf_name(sample_name),
                results=results,
                label_priority=label_priority,
                paths=paths,
                online_status=online_status,
                low_depth_threshold=low_depth_threshold,
                run_id=run_id,
                provenance=provenance,
            )

    window.close()
    sys.exit()