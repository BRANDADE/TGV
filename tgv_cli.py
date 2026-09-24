#!/usr/bin/env python3
import os
import csv
import argparse
import logging
import yaml

from scripts.models.run import Run
from scripts.models.trid import TRID
from scripts.models.dto import AnalysisInput

from scripts.core.vcf_loader import list_vcfs
from scripts.core.vcf_parser import parse_vcf_for_sample
from scripts.core.trid_detector import autodetect_trids
from scripts.core.sequence_utils import reverse_complement
from scripts.core.config_manager import get_safe_config_path
from scripts.core.orchestrator import process_clinical, process_result, process_display

from scripts.bio.clinical_thresholds_loader import load_clinical_thresholds
from scripts.bio.clinical_config_builder import build_clinical_config

# Minimal logging configuration for console output
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

# Optional import for Excel export
try:
    import openpyxl
    HAS_OPENPYXL = True
except ImportError:
    HAS_OPENPYXL = False


# Ordre canonique de référence pour l'organisation des colonnes dans le fichier de sortie
CANONICAL_HEADERS = [
    "trgt_version",
    "bed_version",
    "run_id",
    "sample_id",
    "Locus",
    "Depth",
    "Genotype",
    "Classification",
    "Comments"
]


# ---------------------------------------------------------------------------
# Low Coverage Detection (based on warning symbol)
# ---------------------------------------------------------------------------
def check_low_coverage_flag(display_row):
    """
    Detects the presence of the warning symbol ⚠️ in formatted display row
    fields to flag insufficient coverage.
    """
    fields_to_check = [
        display_row.locus,
        display_row.depth,
        display_row.genotype,
        display_row.classification
    ]
    for field in fields_to_check:
        field_str = str(field)
        if "⚠️" in field_str or "\u26a0" in field_str:
            return True
    return False


# ---------------------------------------------------------------------------
# Load Ataxia panel configuration from config_manager
# ---------------------------------------------------------------------------
def load_panel_ataxie():
    path = get_safe_config_path("buttons_panel.yaml")

    if not os.path.exists(path):
        raise FileNotFoundError(f"Configuration file 'buttons_panel.yaml' not found at: {path}")

    try:
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
    except Exception as e:
        raise RuntimeError(f"Error reading buttons_panel.yaml: {e}")

    if not isinstance(data, dict) or "Ataxie" not in data:
        raise KeyError("The panel section 'Ataxie' is missing from buttons_panel.yaml.")

    return data["Ataxie"]


# ---------------------------------------------------------------------------
# Header Alignment Utilities
# ---------------------------------------------------------------------------
def get_union_headers(existing_headers, new_headers):
    """
    Calcule l'union ordonnée des en-têtes existants et des nouveaux en-têtes
    en se basant sur l'ordre canonique défini.
    """
    union_set = set(existing_headers).union(set(new_headers))
    ordered = [h for h in CANONICAL_HEADERS if h in union_set]
    # Ajout d'éventuelles colonnes inattendues non définies dans CANONICAL_HEADERS
    for h in union_set:
        if h not in ordered:
            ordered.append(h)
    return ordered


# ---------------------------------------------------------------------------
# Flat Files Writers (CSV/TSV) with Header Re-alignment
# ---------------------------------------------------------------------------
def write_csv_tsv(out_path, rows, delimiter):
    file_exists = os.path.exists(out_path) and os.path.getsize(out_path) > 0
    new_headers = list(rows[0].keys())

    if not file_exists:
        # Création d'un nouveau fichier avec uniquement les colonnes fournies
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=new_headers, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(rows)
        return

    # Lecture des en-têtes et données existantes
    existing_rows = []
    existing_headers = []
    with open(out_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        existing_headers = reader.fieldnames if reader.fieldnames else []
        for r in reader:
            existing_rows.append(r)

    union_headers = get_union_headers(existing_headers, new_headers)

    if union_headers == existing_headers:
        # Les colonnes correspondent, nous pouvons simplement ajouter à la suite
        with open(out_path, "a", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=union_headers, delimiter=delimiter, extrasaction='ignore')
            for row in rows:
                formatted_row = {h: row.get(h, "") for h in union_headers}
                writer.writerow(formatted_row)
    else:
        # Reconstitution globale du fichier suite à une modification de la structure des colonnes
        all_rows = []
        for r in existing_rows:
            all_rows.append({h: r.get(h, "") for h in union_headers})
        for r in rows:
            all_rows.append({h: r.get(h, "") for h in union_headers})

        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=union_headers, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(all_rows)


# ---------------------------------------------------------------------------
# Excel File Writer (XLSX) with Header Re-alignment
# ---------------------------------------------------------------------------
def write_xlsx(out_path, rows):
    if not HAS_OPENPYXL:
        raise ImportError(
            "The 'openpyxl' package is required to write Excel files (.xlsx). "
            "Please install it using 'pip install openpyxl' or choose '.csv' / '.tsv' output instead."
        )

    file_exists = os.path.exists(out_path) and os.path.getsize(out_path) > 0
    new_headers = list(rows[0].keys())

    if not file_exists:
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = "TRGT Export"
        ws.append(new_headers)
        for row in rows:
            ws.append([row.get(h, "") for h in new_headers])
        wb.save(out_path)
        return

    wb = openpyxl.load_workbook(out_path)
    ws = wb.active

    existing_headers = []
    if ws.max_row >= 1:
        existing_headers = [ws.cell(row=1, column=col).value for col in range(1, ws.max_column + 1)]
        existing_headers = [h for h in existing_headers if h is not None]

    union_headers = get_union_headers(existing_headers, new_headers)

    if union_headers == existing_headers:
        # Ajout simple des nouvelles lignes dans les colonnes existantes
        for row in rows:
            ws.append([row.get(h, "") for h in union_headers])
    else:
        # Lecture des données existantes pour reconstruction
        existing_rows = []
        for r_idx in range(2, ws.max_row + 1):
            row_val = {}
            for c_idx, h in enumerate(existing_headers):
                val = ws.cell(row=r_idx, column=c_idx + 1).value
                row_val[h] = val if val is not None else ""
            existing_rows.append(row_val)

        # Nettoyage de la feuille de calcul
        ws.delete_rows(1, ws.max_row)

        # Réécriture de l'en-tête unifié
        ws.append(union_headers)

        # Réécriture des anciennes lignes alignées sur le nouvel en-tête
        for r in existing_rows:
            ws.append([r.get(h, "") for h in union_headers])

        # Ajout des nouvelles lignes
        for r in rows:
            ws.append([r.get(h, "") for h in union_headers])

    wb.save(out_path)


# ---------------------------------------------------------------------------
# Write Rows to Output File (Appends/Aggregates to TSV, CSV, or XLSX)
# ---------------------------------------------------------------------------
def write_output_file(out_path, rows):
    if not rows:
        return

    _, ext = os.path.splitext(out_path.lower())

    if ext == ".tsv":
        write_csv_tsv(out_path, rows, delimiter="\t")
    elif ext == ".csv":
        write_csv_tsv(out_path, rows, delimiter=",")
    elif ext == ".xlsx":
        write_xlsx(out_path, rows)
    else:
        raise ValueError(
            f"Unsupported file extension '{ext}'. "
            "Supported formats are '.csv', '.tsv', and '.xlsx'."
        )


# ---------------------------------------------------------------------------
# Export Pipeline for the Entire Run
# ---------------------------------------------------------------------------
def export_run(zip_path, trgt_version, bed_version, run_id, out_tsv, all_loci=False):
    logging.info(f"Opening TRGT archive: {zip_path}")

    vcfs = list_vcfs(zip_path)
    if not vcfs:
        raise FileNotFoundError("The ZIP archive does not contain any .trgt.vcf files.")

    run_name = os.path.basename(zip_path)
    run = Run(name=run_name, vcf_zip=zip_path)

    # Autodetect loci and extract static metadata
    trids, _, _, static_info = autodetect_trids(zip_path)

    # Load clinical thresholds configuration
    thresholds_data = load_clinical_thresholds()
    label_priority = thresholds_data.get("label_priority", {})
    low_depth_threshold = thresholds_data.get("low_depth_threshold", None)

    if not label_priority:
        raise ValueError("The 'label_priority' section is required in clinical_thresholds.yaml.")

    # Configure TRID data structures
    for trid_id in trids:
        t = TRID(trid_id)

        if trid_id in static_info:
            info = static_info[trid_id]
            t.chrom = info["chrom"]
            t.start = info["start"]
            t.end = info["end"]
            t.motifs = info["motifs"]

        trid_yaml_block = thresholds_data.get(trid_id)
        if trid_yaml_block:
            t.clinical = build_clinical_config(trid_id, thresholds_data)

            if t.clinical and t.clinical.orientation.lower() == "rc":
                t.motifs_rc = [reverse_complement(m) for m in t.motifs]

        run.trids[trid_id] = t

    # Set up selected loci list depending on the choice of panel vs. all loci
    if not all_loci:
        panel_ataxie = load_panel_ataxie()
        logging.info(f"Ataxia clinical panel loaded ({len(panel_ataxie)} target loci identified).")
    else:
        logging.info("Global export mode active: all detected loci will be processed.")

    rows = []

    # Sequential processing per sample (patient)
    for vcf_filename in vcfs:
        sample_name = os.path.basename(vcf_filename)
        logging.info(f"Analyzing patient: {sample_name}")

        # Reset sample associations to avoid data collision across iterations
        for trid_id in run.trids:
            run.trids[trid_id].samples = {}

        # Parse the VCF file
        sample_trids = parse_vcf_for_sample(
            zip_path=zip_path,
            vcf_filename=vcf_filename,
            global_trids=run.trids,
        )

        # Map parsed sample objects back to run TRIDs
        for trid_id, sample_obj in sample_trids.items():
            if trid_id in run.trids:
                run.trids[trid_id].samples[vcf_filename] = sample_obj

        # Filter and select loci based on the requested output mode
        trids_present = set(sample_trids.keys())
        if all_loci:
            selected_trids = [t for t in trids if t in trids_present]
        else:
            selected_trids = [t for t in panel_ataxie if t in trids_present]

        if not selected_trids:
            logging.warning(f"No targeted loci identified for sample {sample_name}")
            continue

        trids_global = {trid_id: run.trids[trid_id] for trid_id in selected_trids}
        samples = {trid_id: sample_trids[trid_id] for trid_id in selected_trids}

        paths = {
            "vcf": zip_path,
            "repeat_reads": None,
            "spanning_bam": None,
            "motifs_allele": None,
            "motifs_waterfall": None,
            "meth_allele": None,
            "meth_waterfall": None,
            "genome_fasta": None
        }

        analysis_input = AnalysisInput(
            sample_name=vcf_filename,
            trids=trids_global,
            samples=samples,
            ordered_trids=selected_trids,
            paths=paths,
            label_priority=label_priority
        )

        # Run clinical evaluation and format results
        discordances = process_clinical(analysis_input)
        if discordances:
            logging.warning(f"Motif discordances (BED vs clinical YAML) for {sample_name}: {discordances}")
        process_result(analysis_input)

        # Structure results into TSV format
        for trid_id in selected_trids:
            sample_obj = analysis_input.samples.get(trid_id)
            if not sample_obj or not getattr(sample_obj, "result", None):
                continue

            clinical_cfg = analysis_input.trids[trid_id].clinical
            process_display(sample_obj.result, clinical_cfg, low_depth_threshold)

            r = sample_obj.result
            dr = r.display_row

            # Check for low coverage flags
            is_low = check_low_coverage_flag(dr)
            comments = "Low coverage" if is_low else ""

            # Sanitize sample name for presentation
            clean_sample_id = sample_name.replace(".trgt.vcf", "").replace(".vcf", "")

            # Reconstruction de la ligne de données en ignorant les arguments non spécifiés (None)
            row = {}
            if trgt_version is not None:
                row["trgt_version"] = trgt_version
            if bed_version is not None:
                row["bed_version"] = bed_version
            if run_id is not None:
                row["run_id"] = run_id

            row.update({
                "sample_id": clean_sample_id,
                "Locus": dr.locus,
                "Depth": dr.depth,
                "Genotype": dr.genotype,
                "Classification": dr.classification,
                "Comments": comments,
            })
            rows.append(row)

    if rows:
        write_output_file(out_tsv, rows)
        logging.info(f"Results successfully exported to: {out_tsv}")
    else:
        logging.warning("No rows could be generated for export.")


# ---------------------------------------------------------------------------
# Main Execution Block
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Structured export for TRGT run analysis (command-line utility).")
    parser.add_argument("--zip", required=True, help="Path to the trgt_vcfs.zip archive")
    parser.add_argument("--trgt-version", default=None, help="TRGT version used")
    parser.add_argument("--bed-version", default=None, help="Clinical BED file version used")
    parser.add_argument("--run-id", default=None, help="Unique identifier for the run")
    parser.add_argument("--out", default="tgv_export_ataxie.tsv", help="Output filename (must end in .csv, .tsv, or .xlsx)")
    parser.add_argument("--all-loci", action="store_true", help="Export all detected loci instead of filtering for the clinical (Ataxia) panel only")
    args = parser.parse_args()

    # Pre-validation of output format and dependencies
    _, ext = os.path.splitext(args.out.lower())
    if ext not in [".csv", ".tsv", ".xlsx"]:
        parser.error(
            f"Unsupported output file format '{ext}'. "
            "Please specify a filename ending in '.csv', '.tsv', or '.xlsx'."
        )

    if ext == ".xlsx" and not HAS_OPENPYXL:
        parser.error(
            "The 'openpyxl' package is required for Excel exports (.xlsx) but is not installed. "
            "Please install it (pip install openpyxl) or use '.csv' / '.tsv' output instead."
        )

    try:
        export_run(
            zip_path=args.zip, 
            trgt_version=args.trgt_version, 
            bed_version=args.bed_version, 
            run_id=args.run_id, 
            out_tsv=args.out,
            all_loci=args.all_loci
        )
    except Exception as e:
        logging.error(f"Data export failed: {e}", exc_info=True)


if __name__ == "__main__":
    main()