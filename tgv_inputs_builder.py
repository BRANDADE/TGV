import os
import sys
import argparse
import shutil
import logging
import gzip
import json
import csv
import re
from datetime import datetime
import concurrent.futures
import subprocess
import zipfile
import glob

# Chargement optionnel des librairies
try:
    import json5
except ImportError:
    json5 = None

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# Crée le dossier logs à la racine du script
LOG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
os.makedirs(LOG_DIR, exist_ok=True)

# Crée un fichier de log horodaté
timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
log_file = os.path.join(LOG_DIR, f"tgv_input_builder.{timestamp}.log")


def setup_logging(log_file=None):
    log_format = "%(asctime)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s"

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)

    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)

    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[console_handler, file_handler]
    )


def log_file_only(message):
    logger = logging.getLogger()

    file_handler = None
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            file_handler = handler
            break

    if file_handler is None:
        logger.error("No FileHandler found in logger handlers.")
        return

    record = logger.makeRecord(
        name=logger.name,
        level=logging.INFO,
        fn="",
        lno=0,
        msg=message,
        args=None,
        exc_info=None
    )

    file_handler.handle(record)


def check_fasta_index(fasta, samtools_exe="samtools"):
    fai = os.path.abspath(fasta) + ".fai"
    if not os.path.exists(fai):
        logging.error(f"FASTA index missing: {fai}")
        logging.error(f"Run: {samtools_exe} faidx {fasta}")
        sys.exit(1)


def check_file(file):
    if not os.path.exists(os.path.abspath(file)):
        logging.error(f"File not found: {file}")
        sys.exit(1)


def check_executable(exe, name):
    if os.path.exists(exe):
        return
    if shutil.which(exe) is None:
        logging.error(f"Executable not found: {exe}")
        logging.error(f"Ensure {name} is installed and available in PATH (or specify its path via the options).")
        sys.exit(1)


def check_indexed_file_exists(file_path, index_exts):
    """Vérifie si un fichier et au moins l'un de ses index associés existent et ne sont pas vides."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return False
    for ext in index_exts:
        idx_path_1 = file_path + ext
        if os.path.exists(idx_path_1) and os.path.getsize(idx_path_1) > 0:
            return True
        if file_path.endswith(".bam") and ext == ".bai":
            idx_path_2 = file_path[:-4] + ".bai"
            if os.path.exists(idx_path_2) and os.path.getsize(idx_path_2) > 0:
                return True
    return False


def parse_bed_file(bed_path):
    """Parcourt le fichier BED une seule fois et extrait les IDs uniques et leurs intervalles associés."""
    unique_ids = []
    intervals = []
    seen = set()
    if not os.path.exists(bed_path):
        return unique_ids, intervals

    with open(bed_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, col3 = parts[0], parts[1], parts[2], parts[3]

            r_id = col3
            if "ID=" in col3:
                for sp in col3.split(";"):
                    if sp.startswith("ID="):
                        r_id = sp.split("=")[1]
                        break

            if r_id not in seen:
                seen.add(r_id)
                unique_ids.append(r_id)
                intervals.append((chrom, int(start), int(end), r_id))

    return unique_ids, intervals


def parse_list_samples(tsv_path, default_karyotype="XX"):
    """Lit le TSV (2 ou 3 colonnes) et associe à chaque patient son BAM et son caryotype."""
    samples = {}

    with open(tsv_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 2 or len(parts) > 3:
                logging.error(f"Invalid line in {tsv_path}: {line}")
                logging.error("Expected format: <sample_id> <tab> <bam_path> [<tab> <karyotype>]")
                continue

            sample_id = parts[0]
            bam_path = parts[1]

            karyotype = parts[2] if len(parts) == 3 else default_karyotype
            karyotype = karyotype.upper()
            if karyotype not in ["XX", "XY"]:
                if not os.path.exists(karyotype):
                    logging.warning(f"Invalid karyotype '{karyotype}' for sample '{sample_id}'. Defaulting to '{default_karyotype}'.")
                    karyotype = default_karyotype

            if not os.path.exists(bam_path):
                logging.warning(f"BAM not found for sample '{sample_id}': {bam_path}")
                logging.warning(f"Sample '{sample_id}' will be skipped.")
                continue

            samples[sample_id] = {
                "bam_path": bam_path,
                "karyotype": karyotype
            }

    logging.info(f"{len(samples)} valid samples loaded from {tsv_path}")
    return samples


# --- VALIDATEURS DE PARAMÈTRES TRGT ---

def validate_karyotype(val):
    normalized = val.upper()
    if normalized in ["XX", "XY"]:
        return True, normalized
    if os.path.exists(val):
        return True, val
    return False, "Must be 'XX', 'XY' or a valid path to a file."


def validate_preset(val):
    normalized = val.lower()
    if normalized in ["wgs", "targeted"]:
        return True, normalized
    return False, "Must be 'wgs' or 'targeted'."


def validate_color(val):
    normalized = val.lower()
    if normalized in ["always", "auto", "never"]:
        return True, normalized
    return False, "Must be 'always', 'auto' or 'never'."


def validate_genotyper(val):
    normalized = val.lower()
    if normalized in ["size", "cluster"]:
        return True, normalized
    return False, "Must be 'size' or 'cluster'."


def validate_int(val):
    try:
        ival = int(val)
        if ival >= 0:
            return True, ival
        return False, "Must be a non-negative integer."
    except ValueError:
        return False, "Must be a valid integer."


def validate_bool(val):
    normalized = str(val).lower().strip()
    if normalized in ["true", "yes", "y", "1"]:
        return True, True
    if normalized in ["false", "no", "n", "0"]:
        return True, False
    return False, "Must be a boolean (true/false or yes/no)."


def validate_plot_mode(val):
    normalized = str(val).lower().strip()
    if normalized in ["all", "meth", "allele", "waterfall", "motifs"]:
        return True, normalized
    return False, "Must be 'all', 'meth', 'allele', 'waterfall', or 'motifs'."


def validate_string(val):
    return True, str(val).strip()


VALIDATORS = {
    "karyotype": validate_karyotype,
    "preset": validate_preset,
    "verbose": validate_int,
    "color": validate_color,
    "genotyper": validate_genotyper,
    "flank-len": validate_int,
    "max-depth": validate_int,
    "output-flank-len": validate_int,
}

PARAM_HELP = {
    "karyotype": "Possible values: 'XX', 'XY', or a valid path to a file.",
    "preset": "Possible values: 'wgs', 'targeted'.",
    "verbose": "Possible values:\n          0 = Standard mode (no -v)\n          1 = Verbose mode (-v)\n          2 = Very verbose mode (-vv)",
    "color": "Possible values: 'always', 'auto', 'never'.",
    "genotyper": "Possible values: 'size', 'cluster'.",
    "flank-len": "Possible values: Any positive integer (e.g., 250).",
    "max-depth": "Possible values: Any positive integer (e.g., 250).",
    "output-flank-len": "Possible values: Any positive integer (e.g., 50).",
}

PLOT_VALIDATORS = {
    "plot_mode": validate_plot_mode,
    "verbose": validate_int,
    "squished": validate_bool,
    "font-family": validate_string,
    "color": validate_color,
    "flank-len": validate_int,
    "max-allele-reads": validate_int,
}

PLOT_PARAM_HELP = {
    "plot_mode": "Possible values: 'all', 'meth', 'allele', 'waterfall', or 'motifs'.",
    "verbose": "Possible values:\n          0 = Standard mode (no -v)\n          1 = Verbose mode (-v)\n          2 = Very verbose mode (-vv)",
    "squished": "Possible values: 'true', 'false', 'yes', 'no'.",
    "font-family": "Possible values: Any font family name (e.g., 'Roboto Mono').",
    "color": "Possible values: 'always', 'auto', 'never'.",
    "flank-len": "Possible values: Any positive integer (e.g., 50).",
    "max-allele-reads": "Possible values: Any positive integer (e.g., 100).",
}


def ask_param_validated(key, current_val):
    validator = VALIDATORS.get(key)
    if not validator:
        print(f"  {key} [{current_val}]: ", end="")
        new_val = input().strip()
        return new_val if new_val else current_val

    while True:
        print(f"  {key} [{current_val}]: ", end="")
        log_file_only(f"Prompt (genotype): {key} [{current_val}]")
        new_val = input().strip()
        log_file_only(f"User input (genotype): {new_val}")

        if not new_val:
            return current_val

        is_valid, validated_val = validator(new_val)
        if is_valid:
            log_file_only(f"Value accepted for {key}: {validated_val}")
            return validated_val

        print(f"  [ERROR] Invalid value: {validated_val}")
        print(f"  [HELP]  {PARAM_HELP.get(key, '')}\n")


def ask_plot_param_validated(key, current_val):
    validator = PLOT_VALIDATORS.get(key)
    if not validator:
        print(f"  {key} [{current_val}]: ", end="")
        new_val = input().strip()
        return new_val if new_val else current_val

    while True:
        print(f"  {key} [{current_val}]: ", end="")
        log_file_only(f"Prompt (plot): {key} [{current_val}]")
        new_val = input().strip()
        log_file_only(f"User input (plot): {new_val}")

        if not new_val:
            return current_val

        is_valid, validated_val = validator(new_val)
        if is_valid:
            log_file_only(f"Value accepted for plot {key}: {validated_val}")
            return validated_val

        print(f"  [ERROR] Invalid value: {validated_val}")
        print(f"  [HELP]  {PLOT_PARAM_HELP.get(key, '')}\n")


def build_trgt_command(sample_id, bam_path, karyotype, threads_for_this_sample, args, trgt_params, output_root):
    sample_out = os.path.join(output_root, sample_id)
    os.makedirs(sample_out, exist_ok=True)

    output_prefix = os.path.join(sample_out, f"{sample_id}.trgt")

    cmd = [
        args.trgt,
        "genotype",
        "--genome", args.reference,
        "--reads", bam_path,
        "--repeats", args.bed,
        "--output-prefix", output_prefix,
        "--threads", str(threads_for_this_sample),
        "--karyotype", karyotype
    ]

    EXCLUDED_PARAMS = {"threads", "sample_name", "disable_bam_output", "karyotype"}

    for key, value in trgt_params.items():
        if key in EXCLUDED_PARAMS:
            continue
        if key == "verbose":
            if isinstance(value, int) and value > 0:
                cmd.append("-" + "v" * value)
            continue
        if value is None:
            continue
        if isinstance(value, bool):
            if value:
                cmd.append(f"--{key}")
            continue
        cmd.extend([f"--{key}", str(value)])

    log_file_only(f"TRGT command for {sample_id}: {' '.join(cmd)}")
    return cmd, sample_out


def run_plot_command(plot_cmd, sample_id, repeat_id, plot_type, show):
    """Exécute individuellement une commande de trgt plot et capture ses erreurs."""
    try:
        subprocess.run(plot_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    except subprocess.CalledProcessError as e:
        log_file_only(
            f"Warning: Failed to generate {plot_type}/{show} plot for repeat ID '{repeat_id}' "
            f"(Sample: {sample_id}). Stderr: {e.stderr.strip()}"
        )


def get_plot_combinations(plot_mode="all"):
    normalized = str(plot_mode).lower().strip()

    all_combos = [
        ("allele", "motifs"),
        ("allele", "meth"),
        ("waterfall", "motifs"),
        ("waterfall", "meth")
    ]

    if normalized == "all":
        return all_combos
    elif normalized == "meth":
        return [c for c in all_combos if c[1] == "meth"]
    elif normalized == "motifs":
        return [c for c in all_combos if c[1] == "motifs"]
    elif normalized == "allele":
        return [c for c in all_combos if c[0] == "allele"]
    elif normalized == "waterfall":
        return [c for c in all_combos if c[0] == "waterfall"]
    else:
        return all_combos


def create_sample_zip(directory, zip_path):
    svg_files = glob.glob(os.path.join(directory, "*.trvz.svg"))
    if not svg_files:
        return False
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
        for f in svg_files:
            zf.write(f, os.path.basename(f))
    return True


def create_global_aggregations(output_root, run_name, samples, skip_plots=False):
    logging.info("Starting final global aggregations...")

    if not skip_plots:
        categories = ["motifs_allele", "motifs_waterfall", "meth_allele", "meth_waterfall"]
        for cat in categories:
            zip_name = f"{run_name}-trgt_{cat}.zip"
            zip_path = os.path.join(output_root, zip_name)

            search_pattern = os.path.join(output_root, "*", f"*_{cat}.trvz_alleles.zip")
            sample_zips = glob.glob(search_pattern)

            if sample_zips:
                logging.info(f"Aggregating {len(sample_zips)} zip files into {zip_name}...")
                with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as main_zip:
                    for szip in sample_zips:
                        main_zip.write(szip, os.path.basename(szip))
            else:
                log_file_only(f"No sample-level zips found for category: {cat}")

    bam_zip_name = f"{run_name}-spanning_BAM.zip"
    bam_zip_path = os.path.join(output_root, bam_zip_name)
    bams = glob.glob(os.path.join(output_root, "*", "*.trgt.spanning.sorted.bam"))
    bais = glob.glob(os.path.join(output_root, "*", "*.trgt.spanning.sorted.bam.bai"))
    alt_bais = glob.glob(os.path.join(output_root, "*", "*.trgt.spanning.sorted.bai"))
    all_bam_files = bams + bais + alt_bais

    if all_bam_files:
        logging.info(f"Aggregating {len(all_bam_files)} BAM/BAI files into {bam_zip_name} (Fast Copy mode)...")
        with zipfile.ZipFile(bam_zip_path, 'w', zipfile.ZIP_STORED) as bam_zip:
            for bfile in all_bam_files:
                bam_zip.write(bfile, os.path.basename(bfile))
    else:
        logging.warning("No spanning BAM/BAI files found for aggregation.")

    vcf_zip_name = f"{run_name}-trgt_vcfs.zip"
    vcf_zip_path = os.path.join(output_root, vcf_zip_name)
    vcfs = glob.glob(os.path.join(output_root, "*", "*.trgt.sorted.vcf.gz"))

    if vcfs:
        logging.info(f"Decompressing and aggregating {len(vcfs)} VCF files into {vcf_zip_name}...")
        with zipfile.ZipFile(vcf_zip_path, 'w', zipfile.ZIP_DEFLATED) as vcf_zip:
            for vgz_path in vcfs:
                base_name = os.path.basename(vgz_path)
                vcf_name = base_name[:-3] if base_name.endswith(".gz") else base_name

                with gzip.open(vgz_path, 'rb') as f_in:
                    with vcf_zip.open(vcf_name, 'w') as f_out:
                        shutil.copyfileobj(f_in, f_out)
    else:
        logging.warning("No VCF files found for aggregation.")

    input_bam_zip_name = f"{run_name}-repeat_reads.zip"
    input_bam_zip_path = os.path.join(output_root, input_bam_zip_name)

    input_files_to_zip = []
    for sample_id, info in samples.items():
        bpath = info["bam_path"]
        if os.path.exists(bpath):
            input_files_to_zip.append(bpath)
            bai_1 = bpath + ".bai"
            bai_2 = bpath[:-4] + ".bai" if bpath.endswith(".bam") else None
            if os.path.exists(bai_1):
                input_files_to_zip.append(bai_1)
            elif bai_2 and os.path.exists(bai_2):
                input_files_to_zip.append(bai_2)

    if input_files_to_zip:
        logging.info(f"Aggregating {len(input_files_to_zip)} initial BAM/BAI files into {input_bam_zip_name} (Fast Copy mode)...")
        with zipfile.ZipFile(input_bam_zip_path, 'w', zipfile.ZIP_STORED) as ib_zip:
            for ffile in input_files_to_zip:
                ib_zip.write(ffile, os.path.basename(ffile))
    else:
        logging.warning("No initial BAM/BAI files found for repeat_reads aggregation.")


def run_trgt_genotype(sample_id, bam_path, karyotype, threads_for_this_sample, args, trgt_params, output_root):
    """Exécute la phase de génotypage, tri et indexation pour un échantillon."""
    sample_out = os.path.join(output_root, sample_id)
    os.makedirs(sample_out, exist_ok=True)

    sorted_vcf = os.path.join(sample_out, f"{sample_id}.trgt.sorted.vcf.gz")
    sorted_bam = os.path.join(sample_out, f"{sample_id}.trgt.spanning.sorted.bam")

    vcf_ok = check_indexed_file_exists(sorted_vcf, [".tbi", ".csi"])
    bam_ok = check_indexed_file_exists(sorted_bam, [".bai"])

    if vcf_ok and bam_ok:
        logging.info(f"Sample {sample_id} is already genotyped. Skipping genotyping.")
        return sample_id, sorted_vcf, sorted_bam

    logging.info(f"Starting genotype for sample: {sample_id} with {threads_for_this_sample} threads ({karyotype})")

    bai_path = bam_path + ".bai"
    alt_bai_path = bam_path[:-4] + ".bai" if bam_path.endswith(".bam") else None
    has_index = os.path.exists(bai_path) or (alt_bai_path and os.path.exists(alt_bai_path))

    if not has_index:
        logging.warning(f"Input BAM index (.bai) missing for sample {sample_id}. Indexing {bam_path}...")
        try:
            subprocess.run([args.samtools, "index", "-@", str(threads_for_this_sample), bam_path], check=True)
            logging.info(f"Successfully generated index for input BAM of sample {sample_id}.")
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to generate index for input BAM {bam_path} (exit code {e.returncode})")
            raise e

    cmd, sample_out = build_trgt_command(sample_id, bam_path, karyotype, threads_for_this_sample, args, trgt_params, output_root)

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        logging.info(f"TRGT completed successfully for sample {sample_id}")

        if result.stderr:
            log_file_only(f"TRGT stderr output for {sample_id}:\n{result.stderr}")

        raw_vcf = os.path.join(sample_out, f"{sample_id}.trgt.vcf.gz")
        raw_bam = os.path.join(sample_out, f"{sample_id}.trgt.spanning.bam")

        if os.path.exists(raw_vcf):
            logging.info(f"Sorting and indexing VCF for sample {sample_id}...")
            subprocess.run([
                args.bcftools, "sort",
                "-Oz", "-o", sorted_vcf, raw_vcf
            ], check=True)
            subprocess.run([
                args.bcftools, "index",
                "-t", sorted_vcf
            ], check=True)
            if not args.keep_temp:
                os.remove(raw_vcf)
        else:
            logging.error(f"Raw VCF missing for {sample_id}. Skipping VCF sorting.")

        if os.path.exists(raw_bam):
            logging.info(f"Sorting and indexing BAM for sample {sample_id} with {threads_for_this_sample} threads...")
            subprocess.run([
                args.samtools, "sort",
                "-@", str(threads_for_this_sample),
                "-o", sorted_bam, raw_bam
            ], check=True)
            subprocess.run([
                args.samtools, "index",
                "-@", str(threads_for_this_sample),
                sorted_bam
            ], check=True)
            if not args.keep_temp:
                os.remove(raw_bam)
        else:
            logging.error(f"Raw BAM missing for {sample_id}. Skipping BAM sorting.")

    except subprocess.CalledProcessError as e:
        logging.error(f"TRGT or post-processing failed for sample {sample_id} (exit code {e.returncode})")
        logging.error(f"Error details:\n{e.stderr}")
        raise e
    except Exception as e:
        logging.error(f"Unexpected error running genotype for {sample_id}: {e}")
        raise e

    return sample_id, sorted_vcf, sorted_bam


# --- FONCTIONS DE CALCUL ET GÉNÉRATION DE RAPPORT DE QC ---

def get_bam_metrics(bam_path, samtools_exe="samtools"):
    """Extrait les métriques réelles globales du BAM d'entrée de manière optimisée."""
    num_reads = 0
    lengths = []
    qualities = []

    if not os.path.exists(bam_path):
        return num_reads, 0, "Q0"

    try:
        cmd_idx = [samtools_exe, "idxstats", bam_path]
        res = subprocess.run(cmd_idx, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        for line in res.stdout.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 4:
                num_reads += int(parts[2]) + int(parts[3])
    except Exception:
        try:
            cmd_count = [samtools_exe, "view", "-c", bam_path]
            res = subprocess.run(cmd_count, stdout=subprocess.PIPE, text=True, check=True)
            num_reads = int(res.stdout.strip())
        except Exception:
            num_reads = 0

    cmd_view = [samtools_exe, "view", bam_path]
    try:
        process = subprocess.Popen(cmd_view, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        count = 0
        for line in process.stdout:
            parts = line.split("\t")
            if len(parts) > 10:
                seq = parts[9]
                qual = parts[10].strip()
                lengths.append(len(seq))
                if qual and qual != "*":
                    mean_q = sum(ord(c) - 33 for c in qual) / len(qual)
                    qualities.append(mean_q)
                count += 1
                if count >= 10000:
                    break
        process.stdout.close()
        process.terminate()
        process.wait()
    except Exception as e:
        log_file_only(f"Error sampling BAM reads for metrics: {e}")

    if not lengths:
        return num_reads, 0, "Q0"

    lengths.sort()
    qualities.sort()
    median_length = lengths[len(lengths) // 2]
    median_q_val = qualities[len(qualities) // 2] if qualities else 0

    return num_reads, median_length, f"Q{round(median_q_val)}"


def get_flagstat_metrics(bam_path, samtools_exe="samtools"):
    metrics = {"total": 0, "mapped": 0, "duplicates": 0, "pct_mapped": 0.0, "pct_duplicates": 0.0}
    if not os.path.exists(bam_path):
        return metrics

    try:
        cmd = [samtools_exe, "flagstat", bam_path]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        for line in res.stdout.splitlines():
            stripped = line.strip()
            if "in total" in stripped:
                metrics["total"] = int(stripped.split()[0])
            elif "mapped (" in stripped and "primary mapped" not in stripped:
                metrics["mapped"] = int(stripped.split()[0])
            elif stripped.endswith("duplicates") and "primary duplicates" not in stripped:
                metrics["duplicates"] = int(stripped.split()[0])

        if metrics["total"] > 0:
            metrics["pct_mapped"] = 100 * metrics["mapped"] / metrics["total"]
            metrics["pct_duplicates"] = 100 * metrics["duplicates"] / metrics["total"]
    except Exception as e:
        log_file_only(f"Error running samtools flagstat on {bam_path}: {e}")

    return metrics


def get_ontarget_pct(bam_path, bed_path, mapped_reads, samtools_exe="samtools"):
    if not os.path.exists(bam_path) or not os.path.exists(bed_path) or mapped_reads <= 0:
        return 0.0
    try:
        cmd = [samtools_exe, "view", "-c", "-L", bed_path, bam_path]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        on_target = int(res.stdout.strip())
        return 100 * on_target / mapped_reads
    except Exception as e:
        log_file_only(f"Error running samtools view -c -L on {bam_path}: {e}")
        return 0.0


def extract_original_bedcov(bam_path, bed_path, samtools_exe="samtools"):
    coverages = {}
    if not os.path.exists(bam_path) or not os.path.exists(bed_path):
        return coverages, 0.0

    cmd = [samtools_exe, "bedcov", bed_path, bam_path]
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        total_depth = 0
        total_length = 0

        for line in res.stdout.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 5:
                start = int(parts[1])
                end = int(parts[2])
                col3 = parts[3]
                sum_depth = int(parts[4])

                target_id = col3
                if "ID=" in col3:
                    for sp in col3.split(";"):
                        if sp.startswith("ID="):
                            target_id = sp.split("=")[1]
                            break

                length = end - start
                if length > 0:
                    avg_depth = sum_depth / length
                    coverages[target_id] = avg_depth
                    total_depth += sum_depth
                    total_length += length

        mean_coverage = total_depth / total_length if total_length > 0 else 0.0
        return coverages, mean_coverage
    except Exception as e:
        log_file_only(f"Error running samtools bedcov: {e}")
        return {}, 0.0


def parse_per_allele_dp(fmt_dict):
    """
    Retourne la profondeur PAR ALLELE (une valeur par haplotype), au lieu
    d'une profondeur totale. TRGT expose ceci dans le champ FORMAT 'SD'
    (Spanning read Depth per allele), ex: SD=40,2 -> allèle 1 = 40 reads,
    allèle 2 = 2 reads.

    On ne fait JAMAIS la moyenne ou la somme de ces valeurs pour évaluer la
    qualité d'un génotype : un allèle porteur d'une expansion peut être
    couvert par très peu de reads capables de le "spanner" en entier, alors
    que l'autre allèle (normal, plus court) est couvert par beaucoup plus de
    reads. Moyenner masquerait ce déséquilibre, qui est précisément le signal
    à surveiller.

    Fallback sur DP (profondeur totale non décomposable par allèle) si SD
    est absent du VCF -- dans ce cas on ne peut renvoyer qu'une seule valeur,
    ce qui revient à l'ancien comportement pour ce locus uniquement.
    """
    if 'SD' in fmt_dict and fmt_dict['SD'] not in ('.', ''):
        try:
            return [int(v) for v in fmt_dict['SD'].split(',') if v != '.']
        except ValueError:
            pass
    if 'DP' in fmt_dict:
        try:
            return [int(fmt_dict['DP'])]
        except ValueError:
            pass
    return []


def parse_vcf_for_qc(vcf_path, dp_threshold=10):
    """
    Qualifie chaque locus sur son allèle le MOINS couvert (min_dp), et non
    sur la profondeur totale ou moyenne des deux allèles. Voir
    parse_per_allele_dp() pour le raisonnement complet.
    """
    locus_details = {}
    summary = {"n_loci": 0, "pct_pass": 0.0, "mean_min_allele_dp": 0.0, "flagged_loci": []}

    if not os.path.exists(vcf_path):
        return locus_details, summary

    n_pass = 0
    min_dp_values = []
    flagged = []

    try:
        with gzip.open(vcf_path, 'rt', encoding='utf-8') as f_vcf:
            for line in f_vcf:
                if line.startswith('#'):
                    continue
                vcf_parts = line.strip().split('\t')
                if len(vcf_parts) < 10:
                    continue

                vcf_filt = vcf_parts[6]
                info_field = vcf_parts[7]
                fmt_field = vcf_parts[8]
                sample_field = vcf_parts[9]

                locus_id = extract_trid_from_info(info_field)

                f_keys = fmt_field.split(':')
                f_vals = sample_field.split(':')
                fmt_dict = dict(zip(f_keys, f_vals))

                is_pass = (vcf_filt == '.' or vcf_filt == 'PASS')

                dp_per_allele = parse_per_allele_dp(fmt_dict)
                min_dp = min(dp_per_allele) if dp_per_allele else 0

                locus_details[locus_id] = {
                    "dp_per_allele": dp_per_allele,
                    "min_dp": min_dp,
                    "filter": "PASS" if is_pass else vcf_filt,
                }

                min_dp_values.append(min_dp)
                if is_pass:
                    n_pass += 1
                if (not is_pass) or min_dp < dp_threshold:
                    flagged.append({
                        "repeat_id": locus_id,
                        "filter": vcf_filt,
                        "dp_per_allele": dp_per_allele,
                        "min_dp": min_dp
                    })

        n_total = len(locus_details)
        summary = {
            "n_loci": n_total,
            "pct_pass": (100 * n_pass / n_total) if n_total else 0.0,
            "mean_min_allele_dp": (sum(min_dp_values) / len(min_dp_values)) if min_dp_values else 0.0,
            "flagged_loci": flagged
        }
    except Exception as e:
        log_file_only(f"Failed parsing VCF for diagnostic values ({vcf_path}): {e}")

    return locus_details, summary


def generate_qc_boxplots(read_count_matrix, sample_list, repeat_ids, output_root, run_name):
    png_paths = []
    if not HAS_MATPLOTLIB or not repeat_ids or not sample_list:
        return png_paths

    # Graphique 1 : couverture BAM par échantillon
    try:
        data_per_sample = [
            [read_count_matrix[rid][i] for rid in repeat_ids]
            for i in range(len(sample_list))
        ]
        fig_width = max(6, 0.5 * len(sample_list))
        fig, ax = plt.subplots(figsize=(fig_width, 6))
        
        # Gestion de la compatibilité selon la version de matplotlib installée
        try:
            bp = ax.boxplot(data_per_sample, labels=sample_list, showfliers=True, patch_artist=True)
        except TypeError:
            bp = ax.boxplot(data_per_sample, tick_labels=sample_list, showfliers=True, patch_artist=True)
        
        colors = matplotlib.colormaps['Pastel2']
        num_items = len(sample_list)
        for i, patch in enumerate(bp['boxes']):
            patch.set_facecolor(colors(i / max(1, num_items - 1)))
            patch.set_alpha(0.85)
            patch.set_edgecolor('#1a1015')
            
        for whisker in bp['whiskers']:
            whisker.set_color('#1a1015')
            whisker.set_linestyle('--')
            
        for cap in bp['caps']:
            cap.set_color('#1a1015')
            
        for median in bp['medians']:
            median.set_color('#c8336a')
            median.set_linewidth(2)
            
        for flier in bp['fliers']:
            flier.set_marker('o')
            flier.set_markerfacecolor('#9e1f4f')
            flier.set_alpha(0.5)

        ax.set_ylabel("BAM coverage depth (X)", fontsize=10)
        ax.set_xlabel("Sample", fontsize=10)
        ax.set_title(f"{run_name} - BAM coverage distribution per sample", fontsize=11, fontweight='bold')
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        sample_svg = os.path.join(output_root, f"reads_per_sample_boxplot.svg")
        fig.savefig(sample_svg, format='svg')
        plt.close(fig)
        png_paths.append(sample_svg)
    except Exception as e:
        log_file_only(f"Error generating reads_per_sample boxplot: {e}")

    # Graphique 2 : couverture BAM par locus
    try:
        data_per_locus = [read_count_matrix[rid] for rid in repeat_ids]
        fig_width = max(6, 0.3 * len(repeat_ids))
        fig, ax = plt.subplots(figsize=(fig_width, 6))
        
        # Gestion de la compatibilité selon la version de matplotlib installée
        try:
            bp = ax.boxplot(data_per_locus, labels=repeat_ids, showfliers=True, patch_artist=True)
        except TypeError:
            bp = ax.boxplot(data_per_locus, tick_labels=repeat_ids, showfliers=True, patch_artist=True)
        
        colors = matplotlib.colormaps['Pastel2']
        num_items = len(repeat_ids)
        for i, patch in enumerate(bp['boxes']):
            patch.set_facecolor(colors(i / max(1, num_items - 1)))
            patch.set_alpha(0.85)
            patch.set_edgecolor('#1a1015')
            
        for whisker in bp['whiskers']:
            whisker.set_color('#1a1015')
            whisker.set_linestyle('--')
            
        for cap in bp['caps']:
            cap.set_color('#1a1015')
            
        for median in bp['medians']:
            median.set_color('#c8336a')
            median.set_linewidth(2)
            
        for flier in bp['fliers']:
            flier.set_marker('o')
            flier.set_markerfacecolor('#9e1f4f')
            flier.set_alpha(0.5)

        ax.set_ylabel("BAM coverage depth (X)", fontsize=10)
        ax.set_xlabel("Target locus", fontsize=10)
        ax.set_title(f"{run_name} - BAM coverage distribution per locus", fontsize=11, fontweight='bold')
        plt.xticks(rotation=90)
        plt.tight_layout()
        locus_svg = os.path.join(output_root, f"reads_per_locus_boxplot.svg")
        fig.savefig(locus_svg, format='svg')
        plt.close(fig)
        png_paths.append(locus_svg)
    except Exception as e:
        log_file_only(f"Error generating reads_per_locus boxplot: {e}")

    return png_paths


def compute_sample_status(bam_metrics, vcf_summary, min_mean_cov=10.0, min_pct_pass=90.0):
    flags = []

    if bam_metrics["mean_coverage"] < min_mean_cov:
        flags.append(f"low_mean_coverage(<{min_mean_cov}X)")

    if vcf_summary["pct_pass"] < min_pct_pass:
        flags.append(f"low_pct_pass_loci(<{min_pct_pass}%)")

    if vcf_summary["flagged_loci"]:
        flags.append(f"{len(vcf_summary['flagged_loci'])}_flagged_loci")

    if not flags:
        status = "PASS"
    elif bam_metrics["mean_coverage"] < min_mean_cov / 2:
        status = "FAIL"
    else:
        status = "WARN"

    return status, flags


def extract_trid_from_info(info_field):
    match = re.search(r'TRID=([^;]+)', info_field)
    return match.group(1) if match else "UNKNOWN"


def process_sample_qc(sample_id, bam_path, output_root, args, repeat_ids, samtools_exe):
    """Calcule la QC technique d'un échantillon de manière isolée."""
    num_reads, med_len, med_qual = get_bam_metrics(bam_path, samtools_exe)

    flagstat = get_flagstat_metrics(bam_path, samtools_exe)
    pct_ontarget = get_ontarget_pct(bam_path, args.bed, flagstat["mapped"], samtools_exe)

    coverages, mean_cov = extract_original_bedcov(bam_path, args.bed, samtools_exe)

    total_targets = len(repeat_ids) if repeat_ids else 1
    cov_vals = [coverages.get(rid, 0.0) for rid in repeat_ids]

    pct_10x = (sum(1 for v in cov_vals if v >= 10.0) / total_targets) * 100
    pct_20x = (sum(1 for v in cov_vals if v >= 20.0) / total_targets) * 100
    pct_30x = (sum(1 for v in cov_vals if v >= 30.0) / total_targets) * 100
    pct_low = (sum(1 for v in cov_vals if v < 5.0) / total_targets) * 100

    sorted_vcf = os.path.join(output_root, sample_id, f"{sample_id}.trgt.sorted.vcf.gz")
    vcf_details, vcf_summary = parse_vcf_for_qc(sorted_vcf)

    bam_metrics_dict = {
        "num_reads": num_reads,
        "median_read_length": med_len,
        "median_read_quality": med_qual,
        "mean_coverage": mean_cov,
        "pct_targets_10x": pct_10x,
        "pct_targets_20x": pct_20x,
        "pct_targets_30x": pct_30x,
        "pct_targets_low_cov": pct_low,
        "pct_mapped": flagstat["pct_mapped"],
        "pct_duplicates": flagstat["pct_duplicates"],
        "pct_on_target": pct_ontarget
    }

    status, flags = compute_sample_status(bam_metrics_dict, vcf_summary)

    return {
        "sample_id": sample_id,
        "num_reads": num_reads,
        "median_read_length": med_len,
        "median_read_quality": med_qual,
        "coverages": coverages,
        "mean_cov": mean_cov,
        "pct_10x": pct_10x,
        "pct_20x": pct_20x,
        "pct_30x": pct_30x,
        "pct_low": pct_low,
        "pct_ontarget": pct_ontarget,
        "pct_duplicates": flagstat["pct_duplicates"],
        "bam_metrics_dict": bam_metrics_dict,
        "vcf_summary": vcf_summary,
        "vcf_details": vcf_details,
        "status": status,
        "flags": flags
    }


def generate_qc_report(output_root, run_name, samples, repeat_ids, bed_intervals, args, 
                       samtools_exe="samtools", bcftools_exe="bcftools",
                       cov_thresholds=None, cov_colors=None, cov_labels=None):
    """Génère le package ZIP complet contenant l'ensemble des métriques de QC structurées et visuelles."""
    
    # Initialisation des valeurs de secours si aucun argument n'est fourni
    if cov_thresholds is None:
        cov_thresholds = {"high": 50, "medium": 10}
    if cov_colors is None:
        cov_colors = {"high": "#baffc9", "medium": "#ffdfba", "low": "#ffb3ba"}
    if cov_labels is None:
        cov_labels = {
            "high": "Min Allele Coverage >= {high}X",
            "medium": "Min Allele Coverage {medium}X - {high_minus_1}X",
            "low": "Min Allele Coverage < {medium}X"
        }

    logging.info("Starting Quality Control data extraction...")

    sample_summary_data = []
    interpretation_samples = []
    run_metrics_rows = []
    coverage_matrix = {rid: [] for rid in repeat_ids}
    read_count_matrix = {rid: [] for rid in repeat_ids}
    sample_list_ordered = sorted(list(samples.keys()))

    total_bases_run = 0
    total_reads_run = 0
    all_lengths = []
    all_qualities = []

    sum_reads = 0
    sum_lengths = 0
    sum_qualities = 0
    sum_mean_cov = 0
    sum_pct_10x = 0
    sum_pct_20x = 0
    sum_pct_30x = 0
    sum_pct_low = 0
    sum_pct_ontarget = 0
    sum_pct_dup = 0
    num_samples = len(sample_list_ordered)

    total_threads = max(1, getattr(args, "threads", 1))
    max_parallel_samples = max(1, min(num_samples, total_threads))

    sample_results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_parallel_samples) as executor:
        futures = {
            executor.submit(
                process_sample_qc,
                sample_id,
                samples[sample_id]["bam_path"],
                output_root,
                args,
                repeat_ids,
                samtools_exe
            ): sample_id
            for sample_id in sample_list_ordered
        }
        for future in concurrent.futures.as_completed(futures):
            sample_id = futures[future]
            try:
                sample_results[sample_id] = future.result()
            except Exception as e:
                logging.error(f"QC computation failed for sample {sample_id}: {e}")

    vcf_extracted_data = {}

    for sample_id in sample_list_ordered:
        result = sample_results.get(sample_id)
        if result is None:
            continue

        num_reads = result["num_reads"]
        med_len = result["median_read_length"]
        med_qual = result["median_read_quality"]
        mean_cov = result["mean_cov"]
        coverages = result["coverages"]
        pct_10x = result["pct_10x"]
        pct_20x = result["pct_20x"]
        pct_30x = result["pct_30x"]
        pct_low = result["pct_low"]
        pct_ontarget = result["pct_ontarget"]
        vcf_summary = result["vcf_summary"]
        bam_metrics_dict = result["bam_metrics_dict"]
        status = result["status"]
        flags = result["flags"]

        total_reads_run += num_reads
        total_bases_run += (num_reads * med_len)
        all_lengths.append(med_len)
        try:
            all_qualities.append(int(med_qual[1:]))
        except ValueError:
            all_qualities.append(0)

        for rid in repeat_ids:
            coverage_matrix[rid].append(coverages.get(rid, 0.0))

        run_metrics_rows.append({
            "sample_id": sample_id,
            "num_reads": num_reads,
            "total_bases": num_reads * med_len,
            "median_read_length": med_len,
            "median_read_quality": med_qual
        })

        sample_summary_data.append([
            sample_id,
            f"{num_reads:,}".replace(",", "\u202f"),
            f"{med_len:,}".replace(",", "\u202f"),
            med_qual,
            f"{mean_cov:.2f}",
            f"{pct_10x:.2f}%",
            f"{pct_20x:.2f}%",
            f"{pct_30x:.2f}%",
            f"{pct_low:.2f}%",
            f"{pct_ontarget:.2f}%",
            f"{result['pct_duplicates']:.2f}%"
        ])

        sum_reads += num_reads
        sum_lengths += med_len
        sum_qualities += int(med_qual[1:]) if med_qual.startswith("Q") else 0
        sum_mean_cov += mean_cov
        sum_pct_10x += pct_10x
        sum_pct_20x += pct_20x
        sum_pct_30x += pct_30x
        sum_pct_low += pct_low
        sum_pct_ontarget += pct_ontarget
        sum_pct_dup += result["pct_duplicates"]

        interpretation_samples.append({
            "sample_id": sample_id,
            "status": status,
            "flags": flags,
            "bam_qc": bam_metrics_dict,
            "vcf_qc": {
                "n_loci": vcf_summary["n_loci"],
                "pct_pass_loci": vcf_summary["pct_pass"],
                "mean_min_allele_dp_loci": vcf_summary["mean_min_allele_dp"],
                "flagged_loci": vcf_summary["flagged_loci"]
            }
        })

        vcf_extracted_data[sample_id] = result["vcf_details"]

    for rid in repeat_ids:
        for s_idx, sample_id in enumerate(sample_list_ordered):
            bam_cov_val = coverage_matrix[rid][s_idx]
            read_count_matrix[rid].append(bam_cov_val)

    avg_reads = round(sum_reads / num_samples) if num_samples else 0
    avg_length = round(sum_lengths / num_samples) if num_samples else 0
    avg_quality = f"Q{round(sum_qualities / num_samples)}" if num_samples else "Q0"
    avg_mean_cov = sum_mean_cov / num_samples if num_samples else 0.0
    avg_pct_10x = sum_pct_10x / num_samples if num_samples else 0.0
    avg_pct_20x = sum_pct_20x / num_samples if num_samples else 0.0
    avg_pct_30x = sum_pct_30x / num_samples if num_samples else 0.0
    avg_pct_low = sum_pct_low / num_samples if num_samples else 0.0
    avg_pct_ontarget = sum_pct_ontarget / num_samples if num_samples else 0.0
    avg_pct_dup = sum_pct_dup / num_samples if num_samples else 0.0

    sample_summary_data.append([
        "Sample Average",
        f"{avg_reads:,}".replace(",", "\u202f"),
        f"{avg_length:,}".replace(",", "\u202f"),
        avg_quality,
        f"{avg_mean_cov:.2f}",
        f"{avg_pct_10x:.2f}%",
        f"{avg_pct_20x:.2f}%",
        f"{avg_pct_30x:.2f}%",
        f"{avg_pct_low:.2f}%",
        f"{avg_pct_ontarget:.2f}%",
        f"{avg_pct_dup:.2f}%"
    ])

    med_length_global = all_lengths[len(all_lengths)//2] if all_lengths else 0
    med_qual_global = f"Q{round(sum(all_qualities)/len(all_qualities))}" if all_qualities else "Q0"

    global_report_data = {
        "attributes": [
            {"name": "Total Bases", "value": total_bases_run},
            {"name": "Total Reads", "value": total_reads_run},
            {"name": "Median Read Length (bp)", "value": med_length_global},
            {"name": "Median Read Quality", "value": med_qual_global},
            {"name": "Sample Count", "value": len(samples)},
            {"name": "Target Regions", "value": len(repeat_ids)}
        ],
        "_comment": f"Generated at {datetime.now().isoformat()}"
    }

    for row in run_metrics_rows:
        row["pct_of_run_reads"] = (100 * row["num_reads"] / total_reads_run) if total_reads_run else 0.0
        row["pct_of_run_bases"] = (100 * row["total_bases"] / total_bases_run) if total_bases_run else 0.0

    temp_json = os.path.join(output_root, "target_enrichment_report.json")
    temp_summary_csv = os.path.join(output_root, "sample_summary.csv")
    temp_cov_csv = os.path.join(output_root, "target_cov_by_sample.csv")
    temp_reads_csv = os.path.join(output_root, "target_reads_by_sample.csv")
    temp_run_metrics_csv = os.path.join(output_root, f"{run_name}-run_metrics_by_sample.csv")
    temp_interp_json = os.path.join(output_root, "qc_interpretation.json")

    with open(temp_json, "w", encoding="utf-8") as jf:
        json.dump(global_report_data, jf, indent=4)

    with open(temp_summary_csv, "w", newline="", encoding="utf-8") as scf:
        writer = csv.writer(scf)
        writer.writerow([
            "Sample",
            "Number of Reads",
            "Median Mapped Read Length",
            "Median Mapped Read Quality",
            "Mean Target Coverage",
            "Percent of Targets with \u226510-fold Coverage",
            "Percent of Targets with \u226520-fold Coverage",
            "Percent of Targets with \u226530-fold Coverage",
            "Percent of Targets with Low Coverage (<5X)",
            "Percent of On-Target Reads",
            "Percent of Duplicate Reads"
        ])
        writer.writerows(sample_summary_data)

    with open(temp_cov_csv, "w", newline="", encoding="utf-8") as tcf:
        writer = csv.writer(tcf)
        writer.writerow(["targetName"] + sample_list_ordered)
        for rid in repeat_ids:
            row_vals = []
            for s_idx, sample_id in enumerate(sample_list_ordered):
                bam_val = int(round(coverage_matrix[rid][s_idx]))
                dp_per_allele = vcf_extracted_data.get(sample_id, {}).get(rid, {}).get("dp_per_allele", [])
                vcf_val_str = ",".join(str(v) for v in dp_per_allele) if dp_per_allele else "0"
                row_vals.append(f"{vcf_val_str} ({bam_val})")
            writer.writerow([rid] + row_vals)

    with open(temp_reads_csv, "w", newline="", encoding="utf-8") as trf:
        writer = csv.writer(trf)
        writer.writerow(["targetName"] + sample_list_ordered)
        for rid in repeat_ids:
            writer.writerow([rid] + read_count_matrix[rid])

    with open(temp_run_metrics_csv, "w", newline="", encoding="utf-8") as rmf:
        writer = csv.writer(rmf)
        writer.writerow([
            "Sample", "Total Reads", "Total Bases", "Median Read Length (bp)",
            "Median Read Quality", "% of Run Reads", "% of Run Bases"
        ])
        for row in run_metrics_rows:
            writer.writerow([
                row["sample_id"],
                row["num_reads"],
                row["total_bases"],
                row["median_read_length"],
                row["median_read_quality"],
                f"{row['pct_of_run_reads']:.2f}%",
                f"{row['pct_of_run_bases']:.2f}%"
            ])
        writer.writerow([
            "Run Total", total_reads_run, total_bases_run, med_length_global, med_qual_global, "100.00%", "100.00%"
        ])

    boxplot_pngs = generate_qc_boxplots(read_count_matrix, sample_list_ordered, repeat_ids, output_root, run_name)

    advanced_plots_pngs = []
    if HAS_MATPLOTLIB:
        try:
            quality_png_path = os.path.join(output_root, "genotyping_quality_status.svg")
            samples_keys = sorted(list(vcf_extracted_data.keys()))
            loci_keys = sorted(list({locus for s in samples_keys for locus in vcf_extracted_data[s].keys()}))

            if samples_keys and loci_keys:
                grid = [[0 for _ in range(len(samples_keys))] for _ in range(len(loci_keys))]
                
                # Récupération dynamique des seuils
                high_val = cov_thresholds.get("high", 50)
                medium_val = cov_thresholds.get("medium", 10)

                for s_idx, sample in enumerate(samples_keys):
                    for l_idx, locus in enumerate(loci_keys):
                        l_info = vcf_extracted_data[sample].get(locus)
                        if not l_info:
                            grid[l_idx][s_idx] = 0
                            continue
                        
                        min_dp = l_info["min_dp"]

                        if min_dp >= high_val:
                            grid[l_idx][s_idx] = 2  # Couleur haute
                        elif min_dp >= medium_val:
                            grid[l_idx][s_idx] = 1  # Couleur moyenne
                        else:
                            grid[l_idx][s_idx] = 0  # Couleur basse

                fig_width = max(11, len(samples_keys) * 0.75)
                fig_height = max(9, len(loci_keys) * 0.35)
                fig, ax = plt.subplots(figsize=(fig_width, fig_height))
                from matplotlib.colors import ListedColormap

                # Couleurs dynamiques extraites de la configuration
                color_low = cov_colors.get("low", "#ffb3ba")
                color_medium = cov_colors.get("medium", "#ffdfba")
                color_high = cov_colors.get("high", "#baffc9")

                cmap = ListedColormap([color_low, color_medium, color_high])
                ax.imshow(grid, cmap=cmap, aspect='auto', origin='lower', vmin=0, vmax=2)

                # Configuration des ticks principaux
                ax.set_xticks(range(len(samples_keys)))
                ax.set_xticklabels(samples_keys, rotation=45, ha='right', fontsize=12)
                ax.set_yticks(range(len(loci_keys)))
                ax.set_yticklabels(loci_keys, fontsize=11)

                # Ajout de bordures blanches (grid minor) entre chaque cellule
                ax.set_xticks([x - 0.5 for x in range(len(samples_keys) + 1)], minor=True)
                ax.set_yticks([y - 0.5 for y in range(len(loci_keys) + 1)], minor=True)
                ax.grid(which='minor', color='white', linestyle='-', linewidth=2.5)
                ax.tick_params(which='both', bottom=False, left=False, labelbottom=True, labelleft=True)

                # Annotation textuelle dans les cellules (non gras, taille réduite)
                for s_idx in range(len(samples_keys)):
                    for l_idx in range(len(loci_keys)):
                        l_info = vcf_extracted_data[samples_keys[s_idx]].get(loci_keys[l_idx])
                        if l_info:
                            grid_val = grid[l_idx][s_idx]
                            if grid_val in (0, 1):
                                val_str = "/".join(str(v) for v in l_info["dp_per_allele"])
                                ax.text(s_idx, l_idx, val_str, ha='center', va='center',
                                        color='#2d2025', fontsize=8, fontweight='normal')

                ax.set_title("Locus Coverage Quality Matrix", fontsize=14, fontweight='bold', pad=18)
                
                # Formatage dynamique des étiquettes de la légende
                high_minus_1 = high_val - 1

                def safe_format(template, **kwargs):
                    try:
                        return template.format(**kwargs)
                    except Exception:
                        res = template
                        for k, v in kwargs.items():
                            res = res.replace(f"{{{k}}}", str(v))
                        return res

                label_high = safe_format(
                    cov_labels.get("high", "Min Allele Coverage >= {high}X"),
                    high=high_val, medium=medium_val, high_minus_1=high_minus_1
                )
                label_medium = safe_format(
                    cov_labels.get("medium", "Min Allele Coverage {medium}X - {high_minus_1}X"),
                    high=high_val, medium=medium_val, high_minus_1=high_minus_1
                )
                label_low = safe_format(
                    cov_labels.get("low", "Min Allele Coverage < {medium}X"),
                    high=high_val, medium=medium_val, high_minus_1=high_minus_1
                )

                patch_green = mpatches.Patch(color=color_high, label=label_high)
                patch_orange = mpatches.Patch(color=color_medium, label=label_medium)
                patch_red = mpatches.Patch(color=color_low, label=label_low)

                ax.legend(
                    handles=[patch_green, patch_orange, patch_red],
                    bbox_to_anchor=(1.02, 1),
                    loc='upper left',
                    title="Coverage Status",
                    fontsize=14,
                    title_fontsize=16,
                    handlelength=2.5,
                    handleheight=1.5,
                    labelspacing=0.6,
                    borderpad=1.2
                )

                plt.tight_layout()
                fig.savefig(quality_png_path, format='svg')
                plt.close(fig)
                advanced_plots_pngs.append(quality_png_path)

        except Exception as e:
            log_file_only(f"Failed to generate genotyping_quality_status diagnostic plot: {e}")

    interpretation_report = {
        "run_name": run_name,
        "generated_at": datetime.now().isoformat(),
        "n_samples": num_samples,
        "n_repeat_targets": len(repeat_ids),
        "run_status": "FAIL" if any(s["status"] == "FAIL" for s in interpretation_samples)
                      else ("WARN" if any(s["status"] == "WARN" for s in interpretation_samples) else "PASS"),
        "samples": interpretation_samples
    }
    with open(temp_interp_json, "w", encoding="utf-8") as ijf:
        json.dump(interpretation_report, ijf, indent=2, ensure_ascii=False)

    qc_zip_name = f"{run_name}-QC.zip"
    qz_zip_path = os.path.join(output_root, qc_zip_name)

    logging.info(f"Compressing QC files into final archive: {qc_zip_name}...")
    with zipfile.ZipFile(qz_zip_path, "w", zipfile.ZIP_DEFLATED) as qz:
        qz.write(temp_json, os.path.basename(temp_json))
        qz.write(temp_summary_csv, os.path.basename(temp_summary_csv))
        qz.write(temp_cov_csv, os.path.basename(temp_cov_csv))
        qz.write(temp_reads_csv, os.path.basename(temp_reads_csv))
        qz.write(temp_run_metrics_csv, os.path.basename(temp_run_metrics_csv))
        qz.write(temp_interp_json, os.path.basename(temp_interp_json))
        for png_path in boxplot_pngs:
            qz.write(png_path, os.path.basename(png_path))
        for png_path in advanced_plots_pngs:
            qz.write(png_path, os.path.basename(png_path))

    temp_files_to_clean = [
        temp_json, temp_summary_csv, temp_cov_csv, temp_reads_csv,
        temp_run_metrics_csv, temp_interp_json
    ] + boxplot_pngs + advanced_plots_pngs
    for temp_file in temp_files_to_clean:
        if os.path.exists(temp_file):
            os.remove(temp_file)

    logging.info("QC standalone report ZIP archive generated successfully.")
    return qz_zip_path


def main():
    setup_logging(log_file)

    logging.info("========================================")
    logging.info("Starting TGV Inputs Builder v1.7.0 (Per-Allele QC Mode)")
    logging.info(f"Python version: {sys.version.split()[0]}")
    logging.info(f"Platform: {sys.platform}")
    logging.info(f"Log file: {os.path.abspath(log_file)}")

    parser = argparse.ArgumentParser(description="TGV Inputs Builder")

    parser.add_argument('--trgt', dest='trgt', default='trgt', help="Path to TRGT executable")
    parser.add_argument('--samtools', dest='samtools', default='samtools', help="Path to samtools executable")
    parser.add_argument('--bcftools', dest='bcftools', default='bcftools', help="Path to bcftools executable")
    parser.add_argument('--reference', '-r', dest='reference', required=True)
    parser.add_argument('--bed', '-b', dest='bed', required=True)
    parser.add_argument('--list_samples', '-l', dest='list_samples', required=True)
    parser.add_argument('--name', '-n', dest='run_name', required=True)
    parser.add_argument('--threads', '-t', dest='threads', default=1, type=int)
    parser.add_argument('--non-interactive', action='store_true')
    parser.add_argument('--keep-temp', action='store_true')
    parser.add_argument('--skip-plots', action='store_true')

    args = parser.parse_args()

    check_executable(args.trgt, "TRGT")
    check_executable(args.samtools, "samtools")
    check_executable(args.bcftools, "bcftools")

    check_file(args.reference)
    check_file(args.bed)
    check_file(args.list_samples)
    check_fasta_index(args.reference, args.samtools)

    params_file = os.path.join("configs", "trgt_params.json5")
    if os.path.exists(params_file) and json5 is not None:
        with open(params_file, "r") as f:
            config = json5.load(f)
            trgt_params = config.get("genotype", {})
            plot_params = config.get("plot", {})
            cov_thresholds = config.get("coverage_thresholds", {})
            cov_colors = config.get("coverage_colors", {})
            cov_labels = config.get("coverage_labels", {})


    else:
        logging.warning("configs/trgt_params.json5 missing or json5 module unavailable. Running with empty parameters.")
        trgt_params = {}
        plot_params = {}

    if not args.non_interactive:
        print("\nCurrent TRGT genotype parameters:")
        for k, v in trgt_params.items():
            print(f"  {k}: {v}")

        print("\nDo you want to modify TRGT genotype parameters ? [y/N]: ", end="")
        custom_geno = input().strip().lower()

        if custom_geno == "y":
            print("\n--- Editing TRGT Genotype Parameters ---")
            for key, value in list(trgt_params.items()):
                trgt_params[key] = ask_param_validated(key, value)
            print()

        if not args.skip_plots:
            print("Current TRGT plot parameters:")
            for k, v in plot_params.items():
                print(f"  {k}: {v}")

            print("\nDo you want to modify TRGT plot parameters ? [y/N]: ", end="")
            custom_plot = input().strip().lower()

            if custom_plot == "y":
                print("\n--- Editing TRGT Plot Parameters ---")
                for key, value in list(plot_params.items()):
                    plot_params[key] = ask_plot_param_validated(key, value)
                print()
    else:
        logging.info("Non-interactive mode active.")

    default_karyotype = trgt_params.get("karyotype", "XX")
    samples = parse_list_samples(args.list_samples, default_karyotype=default_karyotype)
    num_samples = len(samples)

    if num_samples == 0:
        logging.error("No valid samples to process. Exiting.")
        sys.exit(1)

    repeat_ids, bed_intervals = parse_bed_file(args.bed)
    num_repeats = len(repeat_ids)
    logging.info(f"Extracted {num_repeats} unique repeat IDs from BED file.")

    output_root = os.path.join("inputs_" + args.run_name)
    os.makedirs(output_root, exist_ok=True)

    # --- CALCUL INTELLIGENT DU PARALLÉLISME ---
    total_threads = max(1, args.threads)
    target_threads_per_job = 4
    max_parallel_jobs = max(1, min(num_samples, total_threads // target_threads_per_job))
    threads_per_job = max(1, total_threads // max_parallel_jobs)

    logging.info(f"Resource Scheduler Plan:")
    logging.info(f"  - Total thread budget: {total_threads}")
    logging.info(f"  - Parallel genotyping processes: {max_parallel_jobs}")
    logging.info(f"  - Threads per process: {threads_per_job}")

    # =========================================================================
    # PHASE 1 : EXECUTION DU GENOTYPAGE EN PARALLELE
    # =========================================================================
    logging.info("Starting Phase 1: Parallel Genotyping...")
    completed_samples = {}  # {sample_id: (sorted_vcf, sorted_bam)}

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_parallel_jobs) as executor:
        futures = {
            executor.submit(
                run_trgt_genotype,
                sample_id,
                info["bam_path"],
                info["karyotype"],
                threads_per_job,
                args,
                trgt_params,
                output_root
            ): sample_id
            for sample_id, info in samples.items()
        }

        for future in concurrent.futures.as_completed(futures):
            sample_id = futures[future]
            try:
                finished_sample, s_vcf, s_bam = future.result()
                completed_samples[finished_sample] = (s_vcf, s_bam)
                logging.info(f"Finished genotyping sample: {finished_sample}")
            except Exception as e:
                logging.error(f"Genotyping execution failed for sample {sample_id}: {e}")

    # =========================================================================
    # PHASE 2 : GENERATION DE TOUS LES GRAPHES (QUEUE PLATE SANS BRANCHE CONCURRENTE ENFANT)
    # =========================================================================
    if not args.skip_plots and repeat_ids and completed_samples:
        logging.info("Starting Phase 2: Generating TRGT plots...")

        plot_tasks = []
        plot_mode = plot_params.get("plot_mode", "all")
        plot_combinations = get_plot_combinations(plot_mode)

        samples_needing_plots = {}
        for sample_id, (sorted_vcf, sorted_bam) in completed_samples.items():
            sample_out = os.path.join(output_root, sample_id)

            expected_zips = []
            for ptype, show in plot_combinations:
                combo_key = f"{show}_{ptype}"
                zip_filename = f"{sample_id}_{combo_key}.trvz_alleles.zip"
                expected_zips.append(os.path.join(sample_out, zip_filename))

            plots_ok = all(os.path.exists(z) and os.path.getsize(z) > 0 for z in expected_zips)
            if plots_ok:
                logging.info(f"Plots for sample {sample_id} already exist. Skipping plot generation.")
            else:
                samples_needing_plots[sample_id] = (sorted_vcf, sorted_bam)

        if samples_needing_plots:
            for sample_id, (sorted_vcf, sorted_bam) in samples_needing_plots.items():
                sample_out = os.path.join(output_root, sample_id)

                for ptype, show in plot_combinations:
                    combo_key = f"{show}_{ptype}"
                    os.makedirs(os.path.join(sample_out, combo_key), exist_ok=True)

                for repeat_id in repeat_ids:
                    for ptype, show in plot_combinations:
                        combo_key = f"{show}_{ptype}"
                        image_name = f"{repeat_id}.trvz.svg"
                        image_path = os.path.join(sample_out, combo_key, image_name)

                        plot_cmd = [
                            args.trgt, "plot",
                            "--genome", args.reference,
                            "--repeats", args.bed,
                            "--vcf", sorted_vcf,
                            "--spanning-reads", sorted_bam,
                            "--repeat-id", repeat_id,
                            "--plot-type", ptype,
                            "--show", show,
                            "--image", image_path
                        ]

                        verbose_val = plot_params.get("verbose", 0)
                        if isinstance(verbose_val, int) and verbose_val > 0:
                            plot_cmd.append("-" + "v" * verbose_val)

                        squished_val = plot_params.get("squished", False)
                        if isinstance(squished_val, str):
                            _, squished_val = validate_bool(squished_val)
                        if squished_val:
                            plot_cmd.append("--squished")

                        if "font-family" in plot_params and plot_params["font-family"]:
                            plot_cmd.extend(["--font-family", str(plot_params["font-family"])])

                        if "flank-len" in plot_params and plot_params["flank-len"] is not None:
                            plot_cmd.extend(["--flank-len", str(plot_params["flank-len"])])

                        if "max-allele-reads" in plot_params and plot_params["max-allele-reads"] is not None:
                            plot_cmd.extend(["--max-allele-reads", str(plot_params["max-allele-reads"])])

                        if "color" in plot_params and plot_params["color"]:
                            plot_cmd.extend(["--color", str(plot_params["color"])])

                        plot_tasks.append((plot_cmd, sample_id, repeat_id, ptype, show))

            if plot_tasks:
                logging.info(f"Submitting {len(plot_tasks)} unified plot tasks using {total_threads} workers...")
                with concurrent.futures.ThreadPoolExecutor(max_workers=total_threads) as plot_executor:
                    plot_futures = {
                        plot_executor.submit(run_plot_command, cmd, s_id, r_id, ptype, show): (s_id, r_id, ptype, show)
                        for cmd, s_id, r_id, ptype, show in plot_tasks
                    }
                    for future in concurrent.futures.as_completed(plot_futures):
                        s_id, r_id, ptype, show = plot_futures[future]
                        try:
                            future.result()
                        except Exception as e:
                            logging.error(
                                f"Unexpected error while generating {ptype}/{show} plot for "
                                f"repeat ID '{r_id}' (Sample: {s_id}): {e}"
                            )
                logging.info("All plots generated successfully.")

            # =========================================================================
            # PHASE 3 : COMPRESSION INDIVIDUELLE (ZIP PAR PATIENT ET PAR CATEGORIE)
            # =========================================================================
            logging.info("Starting Phase 3: Zipping individual category plot directories...")
            for sample_id in samples_needing_plots.keys():
                sample_out = os.path.join(output_root, sample_id)

                for ptype, show in plot_combinations:
                    combo_key = f"{show}_{ptype}"
                    combo_dir = os.path.join(sample_out, combo_key)
                    zip_filename = f"{sample_id}_{combo_key}.trvz_alleles.zip"
                    zip_filepath = os.path.join(sample_out, zip_filename)

                    if os.path.exists(combo_dir):
                        zip_success = create_sample_zip(combo_dir, zip_filepath)

                        if zip_success:
                            if not args.keep_temp:
                                try:
                                    shutil.rmtree(combo_dir)
                                except Exception as e:
                                    logging.warning(f"Could not clean up directory {combo_dir}: {e}")
                        else:
                            try:
                                shutil.rmtree(combo_dir)
                            except Exception:
                                pass
                            log_file_only(f"No SVG plots found for {sample_id} ({combo_key}). Skipped zipping.")

    # =========================================================================
    # PHASE 3.5 : GENERATION DU RAPPORT DE QC AUTONOME (EXPORT ZIP)
    # =========================================================================
    try:
        generate_qc_report(
            output_root, args.run_name, samples, repeat_ids, bed_intervals, args, 
            samtools_exe=args.samtools, bcftools_exe=args.bcftools,
            cov_thresholds=cov_thresholds, cov_colors=cov_colors, cov_labels=cov_labels
        )
    except Exception as e:
        logging.error(f"Failed to generate Quality Control Report: {e}")

    # =========================================================================
    # PHASE 4 : AGREGATION GLOBALE DES RESULTATS
    # =========================================================================
    create_global_aggregations(output_root, args.run_name, samples, skip_plots=args.skip_plots)

    # =========================================================================
    # PHASE 5 : NETTOYAGE DES DOSSIERS DE SAMPLES INDIVIDUELS
    # =========================================================================
    if not args.keep_temp:
        logging.info("Cleaning up individual sample directories to leave only global archives...")
        for sample_id in samples.keys():
            sample_dir = os.path.join(output_root, sample_id)
            if os.path.exists(sample_dir):
                try:
                    shutil.rmtree(sample_dir)
                    log_file_only(f"Removed individual sample directory: {sample_dir}")
                except Exception as e:
                    logging.warning(f"Could not remove directory {sample_dir}: {e}")
        logging.info("Individual sample directories cleanup completed successfully.")
    else:
        logging.info("Keeping individual sample directories (--keep-temp is active).")

    logging.info("TGV Inputs Builder completed successfully.")


if __name__ == "__main__":
    main()