# scripts/core/i18n.py
import json
import os
import subprocess
import sys

LANG_BUTTON_LABELS = {
    "en": "Language (EN)",
    "fr": "Langue (FR)"
}

TRANSLATIONS = {
    "en": {
        # --- IGV HTML & Bandeau ---
        "Visualisation IGV": "IGV Visualization",
        "MÉTADONNÉES LOCUS": "LOCUS METADATA",
        "Métadonnées Locus": "Locus Metadata",
        "ALLÈLE 1": "ALLELE 1",
        "ALLÈLE 2": "ALLELE 2",
        "Allèle 1": "Allele 1",
        "Allèle 2": "Allele 2",
        "Pureté :": "Purity:",
        "Pureté": "Purity",
        "Génotype :": "Genotype:",
        "Profondeur :": "Depth:",
        "Méthylation :": "Methylation:",
        "Patient :": "Patient:",
        "ID Patient :": "Patient ID:",
        "Locus :": "Locus:",
        "Non classifié": "Unclassified",
        "Non classé": "Unclassified",
        "Motif :": "Motif:",

        # --- HTML Table Export (html_report.py) ---
        "Résultats TGV": "TGV Results",
        "Commentaires": "Comments",
        "Actions": "Actions",
        "Couverture faible": "Low coverage",
        "Ajouter un commentaire...": "Add a comment...",
        "Inverser Allèle 1 / Allèle 2": "Swap Allele 1 / Allele 2",
        "Tout cocher": "Check all",
        "Tout décocher": "Uncheck all",
        "Copier": "Copy",
        "copier": "copy",
        "Copier la sélection": "Copy selection",
        "Copier le tableau": "Copy table",
        "Copier la séquence": "Copy sequence",
        "Copié": "Copied",
        "Copié !": "Copied!",
        "Aucune ligne n'est sélectionnée pour la copie.": "No row selected for copy.",
        "Erreur lors de la copie : ": "Error during copy: ",

        # --- QC HTML Report (make_report.py) ---
        "Target Enrichment Report": "Target Enrichment Report",
        "Run Metrics Summary": "Run Metrics Summary",
        "Run Metrics": "Run Metrics",
        "Global Attributes": "Global Attributes",
        "Copy TSV": "Copy TSV",
        "Copier TSV": "Copy TSV",
        "QC Distribution Plots": "QC Distribution Plots",
        "Distribution Plots": "Distribution Plots",
        "Read Count Distribution Per Sample": "Read Count Distribution Per Sample",
        "Read Count Distribution Per Locus": "Read Count Distribution Per Locus",
        "Locus Coverage Quality Matrix": "Locus Coverage Quality Matrix",
        "Locus Coverage Quality": "Locus Coverage Quality",
        "Sample Summary Dashboard": "Sample Summary Dashboard",
        "Sample Summary": "Sample Summary",
        "Metrics & Genotyping Quality": "Metrics & Genotyping Quality",
        "Global Target Coverage Matrix": "Global Target Coverage Matrix",
        "Target Coverage": "Target Coverage",
        "Save Report": "Save Report",
        "Generated on": "Generated on",
        "Coverage per Sample & Locus (Transposed Matrix - Format: VCF depth per allele (BAM physical coverage))": "Coverage per Sample & Locus (Transposed Matrix - Format: VCF depth per allele (BAM physical coverage))",
        "Interactive Viewport": "Interactive Viewport",
        "Zoom In": "Zoom In",
        "Zoom Out": "Zoom Out",
        "Reset View": "Reset View",
        "Sample ID": "Sample ID",
        "Sample": "Sample",
        "Reads": "Reads",
        "Length (bp)": "Length (bp)",
        "Quality": "Quality",
        "Mean Cov": "Mean Cov",
        "PASS Loci %": "PASS Loci %",
        "Flagged Loci": "Flagged Loci",
        "On-Target %": "On-Target %",
        "Dup %": "Dup %",
        "Impossible d'ouvrir le rapport global :": "Unable to open global report:",
        "Rapport HTML PacBio TRGT & TGV QC": "PacBio TRGT & TGV QC HTML Report",
        "Dossier du run ou fichier ZIP contenant les QC": "Run folder or ZIP archive containing QC files",

        # --- Language Selection ---
        "Sélection de la langue": "Language Selection",
        "Redémarrage de l'application en cours...": "Restarting application, please wait...",

        # --- Window Titles & PySimpleGUI Buttons & Labels ---
        "Résultats pour": "Results for",
        "TGV  TRGT Global Viewer": "TGV  TRGT Global Viewer",
        "Locus": "Locus",
        "Profondeur": "Depth",
        "Taille (bp)": "Size (bp)",
        "Motifs": "Motifs",
        "Génotype": "Genotype",
        "Classification": "Classification",
        "Allèle 1 - Répétition": "Allele 1 - Repeat",
        "Allèle 2 - Répétition": "Allele 2 - Repeat",
        "Auto :": "Auto:",
        "Appliquée :": "Applied:",
        "Appliqué :": "Applied:",
        "Plots :": "Plots:",
        "Détails": "Details",
        "Valider": "Validate",
        "Validate": "Validate",
        "Réinitialiser auto": "Reset auto",
        "Reset auto": "Reset auto",
        "Réinitialiser": "Reset",
        "Séquence allèle 1": "Allele 1 sequence",
        "Séquence allèle 2": "Allele 2 sequence",
        "Ouvrir IGV": "Open IGV",
        "Open IGV": "Open IGV",
        "Ouvrir plot": "Open plot",
        "Open plot": "Open plot",
        "Export Data": "Export Data",
        "Fermer": "Close",
        "Close": "Close",
        "Browse": "Browse",
        "Parcourir": "Browse",
        "Panels": "Panels",
        "Lancer l'analyse": "Launch Analysis",
        "Launch Analysis": "Launch Analysis",
        "Rapport QC": "QC Report",
        "QC Report": "QC Report",
        "Sélection du fichier TRGT (trgt_vcfs.zip)": "Select TRGT file (trgt_vcfs.zip)",
        "Sélection du génome de référence (.fa / .fasta) [Optionnel]": "Select reference genome (.fa / .fasta) [Optional]",
        "Patient à analyser": "Patient to analyze",
        "Recherche locus": "Locus search",
        "Ajouter": "Add",
        "Add": "Add",
        "Locus sélectionnés": "Selected Loci",
        "Retirer": "Remove",
        "Remove": "Remove",
        "Actions sur les locus": "Actions on loci",
        "Inverser les allèles": "Swap alleles",
        "Tout sélectionner": "Select all",
        "Tout désélectionner": "Deselect all",
        "Annuler": "Cancel",
        "Cancel": "Cancel",
        "Quitter": "Exit",
        "Exit": "Exit",
        "OK": "OK",

        # --- Dynamic Text Block ---
        "=== Informations générales ===": "=== General information ===",
        "=== Allèle 1 ===": "=== Allele 1 ===",
        "=== Allèle 2 ===": "=== Allele 2 ===",
        "Motifs TRGT :": "TRGT Motifs:",
        "Génotype auto :": "Auto genotype:",
        "Génotype final :": "Final genotype:",
        "Classification auto :": "Auto classification:",
        "Classification finale :": "Final classification:",
        "Motifs cliniques :": "Clinical motifs:",
        "Répétition :": "Repeat:",
        "Segmentation :": "Segmentation:",
        "Séquence :": "Sequence:",
        "[cliquer bouton]": "[click button]",

        # --- Errors & Status ---
        "Veuillez sélectionner un ZIP et un patient": "Please select a ZIP file and a patient",
        "Veuillez sélectionner un fichier ZIP et un patient": "Please select a ZIP file and a patient",
        "Veuillez sélectionner une archive ZIP et un patient": "Please select a ZIP archive and a patient",
        "Veuillez sélectionner un fichier TRGT valide": "Please select a valid TRGT file",
        "Veuillez sélectionner au moins un locus": "Please select at least one locus",
        "Aucun locus sélectionné": "No locus selected",
        "Aucun patient sélectionné": "No patient selected",
        "Locus non trouvé": "Locus not found",
        "ne sont pas présents dans le VCF": "are not present in the VCF",
        "Certains gènes du panel": "Some panel genes",
        "absent du VCF": "missing from the VCF",
        "absents du VCF": "missing from the VCF",
        "Attention": "Warning",
        "Erreur": "Error",
        "Information": "Information",
        "Succès": "Success",
        "Sain": "Normal",
        "Prémutation": "Premutation",
        "Pathogène": "Pathogenic",
    },
    "fr": {
        "Copy": "Copier",
        "copy": "copier",
        "Copy TSV": "Copier TSV",
        "Copy selection": "Copier la sélection",
        "Copy table": "Copier le tableau",
        "Copy sequence": "Copier la séquence",
        "Copied": "Copié",
        "Copied!": "Copié !",
        "Browse": "Parcourir",
        "FileBrowse": "Parcourir",
        "Launch Analysis": "Lancer l'analyse",
        "QC Report": "Rapport QC",
        "Open IGV": "Ouvrir IGV",
        "Open plot": "Ouvrir plot",
        "Export Data": "Export de données",
        "Close": "Fermer",
        "Validate": "Valider",
        "Reset auto": "Réinitialiser auto",
        "Reset": "Réinitialiser",
        "Cancel": "Annuler",
        "Exit": "Quitter",
        "Select all": "Tout sélectionner",
        "Deselect all": "Tout désélectionner",
        "Target Enrichment Report": "Rapport d'enrichissement de cibles",
        "Run Metrics Summary": "Résumé des métriques du run",
        "Run Metrics": "Métriques du Run",
        "Global Attributes": "Attributs globaux",
        "QC Distribution Plots": "Graphiques de distribution QC",
        "Distribution Plots": "Graphiques de distribution",
        "Read Count Distribution Per Sample": "Distribution des lectures par échantillon",
        "Read Count Distribution Per Locus": "Distribution des lectures par locus",
        "Locus Coverage Quality Matrix": "Matrice de qualité de couverture par locus",
        "Locus Coverage Quality": "Qualité de couverture des locus",
        "Sample Summary Dashboard": "Tableau de bord des échantillons",
        "Sample Summary": "Résumé des échantillons",
        "Metrics & Genotyping Quality": "Métriques & Qualité de génotypage",
        "Global Target Coverage Matrix": "Matrice globale de couverture des cibles",
        "Target Coverage": "Couverture des cibles",
        "Save Report": "Enregistrer le rapport",
        "Generated on": "Généré le",
        "Coverage per Sample & Locus (Transposed Matrix - Format: VCF depth per allele (BAM physical coverage))": "Couverture par échantillon et locus (Matrice transposée - Format : Profondeur VCF par allèle (couverture physique BAM))",
        "Interactive Viewport": "Zone interactive",
        "Zoom In": "Zoom avant",
        "Zoom Out": "Zoom arrière",
        "Reset View": "Réinitialiser la vue",
        "Sample ID": "ID Échantillon",
        "Sample": "Échantillon",
        "Reads": "Lectures",
        "Length (bp)": "Taille (pb)",
        "Quality": "Qualité",
        "Mean Cov": "Couv. Moyenne",
        "PASS Loci %": "% Loci PASS",
        "Flagged Loci": "Loci avec alertes",
        "On-Target %": "% Sur cible",
        "Dup %": "% Doublons",
        "Impossible d'ouvrir le rapport global :": "Impossible d'ouvrir le rapport global :",
        "Rapport HTML PacBio TRGT & TGV QC": "Rapport HTML PacBio TRGT & TGV QC",
        "Dossier du run ou fichier ZIP contenant les QC": "Dossier du run ou fichier ZIP contenant les QC",
        "Sélection de la langue": "Sélection de la langue",
        "Redémarrage de l'application en cours...": "Redémarrage de l'application en cours...",
    }
}

# --- Configuration du chemin vers configs/tgv_lang.json ---
if getattr(sys, 'frozen', False):
    BASE_DIR = os.path.dirname(sys.executable)
else:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    BASE_DIR = os.path.abspath(os.path.join(script_dir, "..", ".."))

CONFIG_DIR = os.path.join(BASE_DIR, "configs")
os.makedirs(CONFIG_DIR, exist_ok=True)
CONFIG_FILE = os.path.join(CONFIG_DIR, "tgv_lang.json")


def get_current_lang():
    if os.path.exists(CONFIG_FILE):
        try:
            with open(CONFIG_FILE, "r", encoding="utf-8") as f:
                return json.load(f).get("language", "en")
        except Exception:
            pass
    return "en"


CURRENT_LANG = get_current_lang()


def set_language(lang_code):
    global CURRENT_LANG
    CURRENT_LANG = lang_code
    try:
        os.makedirs(CONFIG_DIR, exist_ok=True)
        with open(CONFIG_FILE, "w", encoding="utf-8") as f:
            json.dump({"language": lang_code}, f, indent=4)
    except Exception:
        pass


def tr(text: str) -> str:
    if not isinstance(text, str):
        return text

    clean_key = text.strip()
    clean_no_punct = clean_key.rstrip(".:! ")

    lang_dict = TRANSLATIONS.get(CURRENT_LANG, {})

    if clean_key in lang_dict:
        return lang_dict[clean_key]
    
    if clean_no_punct in lang_dict:
        suffix = clean_key[len(clean_no_punct):]
        return lang_dict[clean_no_punct] + suffix

    if CURRENT_LANG == "fr":
        return text

    if clean_key.startswith("Résultats pour"):
        sample_name = clean_key.replace("Résultats pour", "").strip()
        prefix = lang_dict.get("Résultats pour", "Results for")
        return f"{prefix} {sample_name}"

    res = text
    for fr_str, target_str in lang_dict.items():
        if fr_str in res:
            res = res.replace(fr_str, target_str)
    return res


def tr_multiline_block(text: str) -> str:
    if not isinstance(text, str):
        return text
    lang_dict = TRANSLATIONS.get(CURRENT_LANG, {})
    for fr_str, target_str in lang_dict.items():
        if fr_str in text:
            text = text.replace(fr_str, target_str)
    return text


def restart_application(current_window=None):
    """Relance l'application sous Linux ou en mode script Python."""
    if current_window:
        try:
            current_window.close()
        except Exception:
            pass

    if getattr(sys, 'frozen', False):
        cmd = [sys.executable] + sys.argv[1:]
        working_dir = os.path.dirname(sys.executable)
    else:
        cmd = [sys.executable] + sys.argv
        working_dir = BASE_DIR

    # Sous Linux : lancement propre et détaché
    if sys.platform != "win32":
        subprocess.Popen(cmd, cwd=working_dir, start_new_session=True)
    else:
        subprocess.Popen(cmd, cwd=working_dir, shell=True)

    os._exit(0)


# Interceptions PySimpleGUI
try:
    import PySimpleGUI as sg

    def show_language_modal(parent_window=None):
        curr = get_current_lang()
        layout = [
            [sg.Text(tr("Sélection de la langue"), font=("Helvetica", 11, "bold"))],
            [sg.Radio("English (EN)", "LANG_R", default=(curr == "en"), key="r_en")],
            [sg.Radio("Français (FR)", "LANG_R", default=(curr == "fr"), key="r_fr")],
            [sg.Push(), sg.Button("OK", key="-LANG_SAVE-", size=(8, 1)), sg.Button(tr("Fermer"), key="-LANG_CANCEL-", size=(8, 1))]
        ]
        win = sg.Window("Language", layout, modal=True, keep_on_top=True)
        ev, val = win.read(close=True)
        
        if ev == "-LANG_SAVE-":
            new_l = "fr" if val.get("r_fr") else "en"
            set_language(new_l)

            # CAS WINDOWS (Exécutable .exe) : Message clair et fermeture propre
            if sys.platform == "win32" and getattr(sys, 'frozen', False):
                msg = (
                    "Langue modifiée avec succès !\nVeuillez redémarrer l'application pour appliquer les changements."
                    if new_l == "fr" else
                    "Language changed successfully!\nPlease restart the application to apply the changes."
                )
                sg.popup_ok(msg, title="Langue / Language", keep_on_top=True)
                if parent_window:
                    try:
                        parent_window.close()
                    except Exception:
                        pass
                sys.exit(0)
            
            # CAS LINUX OU SCRIPT PYTHON : Redémarrage automatique immédiat
            else:
                restart_application(parent_window)

    # 1. Patch de l'initialisation des boutons
    orig_button_init = sg.Button.__init__
    def patched_button_init(self, *args, **kwargs):
        args_list = list(args)
        orig_text = None
        if len(args_list) > 0 and isinstance(args_list[0], str):
            orig_text = args_list[0]
            if "key" not in kwargs or kwargs["key"] is None:
                kwargs["key"] = orig_text
            args_list[0] = tr(orig_text)
        elif "button_text" in kwargs and isinstance(kwargs["button_text"], str):
            orig_text = kwargs["button_text"]
            if "key" not in kwargs or kwargs["key"] is None:
                kwargs["key"] = orig_text
            kwargs["button_text"] = tr(orig_text)
        orig_button_init(self, *args_list, **kwargs)
    sg.Button.__init__ = patched_button_init

    # 2. Patch dynamique des mises à jour de boutons (.update / .Update)
    def wrap_button_update(orig_fn):
        def patched(self, *args, **kwargs):
            args_list = list(args)
            if len(args_list) > 0 and isinstance(args_list[0], str):
                args_list[0] = tr(args_list[0])
            elif "text" in kwargs and isinstance(kwargs["text"], str):
                kwargs["text"] = tr(kwargs["text"])
            elif "button_text" in kwargs and isinstance(kwargs["button_text"], str):
                kwargs["button_text"] = tr(kwargs["button_text"])
            return orig_fn(self, *args_list, **kwargs)
        return patched

    if hasattr(sg, 'Button'):
        sg.Button.update = wrap_button_update(sg.Button.update)
        if hasattr(sg.Button, 'Update'):
            sg.Button.Update = wrap_button_update(sg.Button.Update)

    # 3. Patch des éléments d'interface textuels
    def patch_simple(target_cls, arg_idx=0, kw="text"):
        orig_init = target_cls.__init__
        def new_init(self, *args, **kwargs):
            args_list = list(args)
            if len(args_list) > arg_idx and isinstance(args_list[arg_idx], str):
                args_list[arg_idx] = tr(args_list[arg_idx])
            elif kw in kwargs and isinstance(kwargs[kw], str):
                kwargs[kw] = tr(kwargs[kw])
            orig_init(self, *args_list, **kwargs)
        target_cls.__init__ = new_init

    patch_simple(sg.Text, 0, "text")
    patch_simple(sg.Frame, 0, "title")
    patch_simple(sg.Tab, 0, "title")
    patch_simple(sg.Checkbox, 0, "text")
    patch_simple(sg.Radio, 0, "text")

    # 4. Patch des mises à jour génériques (.update)
    def wrap_update_method(orig_update_fn):
        def patched_update(self, *args, **kwargs):
            args_list = list(args)
            if len(args_list) > 0 and isinstance(args_list[0], str):
                args_list[0] = tr(args_list[0])
            elif "value" in kwargs and isinstance(kwargs["value"], str):
                kwargs["value"] = tr(kwargs["value"])
            elif "text" in kwargs and isinstance(kwargs["text"], str):
                kwargs["text"] = tr(kwargs["text"])
            return orig_update_fn(self, *args_list, **kwargs)
        return patched_update

    sg.Text.update = wrap_update_method(sg.Text.update)
    sg.Element.update = wrap_update_method(sg.Element.update)
    if hasattr(sg, 'StatusBar'):
        sg.StatusBar.update = wrap_update_method(sg.StatusBar.update)

    orig_table_init = sg.Table.__init__
    def patched_table_init(self, *args, **kwargs):
        if "headings" in kwargs and isinstance(kwargs["headings"], list):
            kwargs["headings"] = [tr(h) for h in kwargs["headings"]]
        elif len(args) > 1 and isinstance(args[1], list):
            args_list = list(args)
            args_list[1] = [tr(h) for h in args_list[1]]
            args = tuple(args_list)
        orig_table_init(self, *args, **kwargs)
    sg.Table.__init__ = patched_table_init

    orig_multiline_update = sg.Multiline.update
    def patched_multiline_update(self, value=None, *args, **kwargs):
        if value is not None and isinstance(value, str):
            value = tr_multiline_block(value)
        return orig_multiline_update(self, value, *args, **kwargs)
    sg.Multiline.update = patched_multiline_update

    def wrap_popup(orig_func):
        def wrapper(*args, **kwargs):
            new_args = [tr(str(a)) if isinstance(a, str) else a for a in args]
            if "title" in kwargs and isinstance(kwargs["title"], str):
                kwargs["title"] = tr(kwargs["title"])
            return orig_func(*new_args, **kwargs)
        return wrapper

    for attr in dir(sg):
        if attr.startswith("popup"):
            try:
                fn = getattr(sg, attr)
                if callable(fn):
                    setattr(sg, attr, wrap_popup(fn))
            except Exception:
                pass

    orig_window_init = sg.Window.__init__
    def patched_window_init(self, title, layout=None, *args, **kwargs):
        title_tr = tr(title)
        is_modal = kwargs.get("modal", False)
        if layout is not None and isinstance(layout, list) and not is_modal and title not in ("Language", "TGV"):
            current_btn_label = LANG_BUTTON_LABELS.get(CURRENT_LANG, "Language (EN)")
            lang_bar = [
                sg.Push(),
                sg.Button(current_btn_label, key="-TGV_OPEN_LANG_MODAL-", size=(14, 1), button_color=('white', '#2c3e50'), pad=(2, 2)),
            ]
            layout = [[lang_bar]] + layout
        orig_window_init(self, title_tr, layout, *args, **kwargs)

    orig_window_read = sg.Window.read
    def patched_window_read(self, *args, **kwargs):
        while True:
            event, values = orig_window_read(self, *args, **kwargs)
            if event == "-TGV_OPEN_LANG_MODAL-":
                show_language_modal(self)
                continue
            return event, values

    sg.Window.__init__ = patched_window_init
    sg.Window.read = patched_window_read

except ImportError:
    pass