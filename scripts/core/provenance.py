"""
Provenance des résultats : versions (TRGT, TGV), catalogue, configuration clinique.
"""
import hashlib
import json
import logging
import os
import shlex
import zipfile

from scripts.core.config_manager import get_safe_config_path
from scripts.core.version import __version__, build_commit

RUN_MANIFEST = "run_manifest.json"


def read_vcf_header(zip_path, vcf_name):
    """
    En-tête TRGT d'un VCF : ##trgtVersion et ##trgtCommand (src/trgt/writers/write_vcf.rs).
    Retourne {"trgt_version", "trgt_command", "catalog"} (chaînes vides si absents).
    """
    info = {"trgt_version": "", "trgt_command": "", "catalog": ""}
    with zipfile.ZipFile(zip_path, "r") as z:
        with z.open(vcf_name) as f:
            for raw in f:
                line = raw.decode("utf-8", errors="replace").rstrip("\n")
                if not line.startswith("##"):
                    break
                if line.startswith("##trgtVersion="):
                    info["trgt_version"] = line.split("=", 1)[1].strip()
                elif line.startswith("##trgtCommand="):
                    info["trgt_command"] = line.split("=", 1)[1].strip()

    info["catalog"] = catalog_from_command(info["trgt_command"])
    return info


def catalog_from_command(command):
    """Nom du catalogue passé à `trgt genotype --repeats/-b`."""
    if not command:
        return ""
    try:
        tokens = shlex.split(command)
    except ValueError:
        tokens = command.split()
    for i, token in enumerate(tokens):
        if token in ("--repeats", "-b") and i + 1 < len(tokens):
            return os.path.basename(tokens[i + 1])
        if token.startswith("--repeats="):
            return os.path.basename(token.split("=", 1)[1])
    return ""


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def thresholds_sha256():
    """SHA-256 du clinical_thresholds.yaml effectivement utilisé."""
    try:
        return file_sha256(get_safe_config_path("clinical_thresholds.yaml"))
    except OSError as e:
        logging.warning(f"Cannot hash clinical_thresholds.yaml: {e}")
        return ""


def format_source(source):
    """Champ 'source' d'un bloc YAML → 'HAS - Volet 1 (02-2025); Benkirane et al. (2025)'."""
    if not isinstance(source, dict):
        return ""
    entries = []
    if "document" in source:
        entries.append((source.get("document"), source.get("version")))
    for doc in source.get("documents") or []:
        if isinstance(doc, dict):
            entries.append((doc.get("name"), doc.get("version")))
    return "; ".join(f"{name} ({version})" if version else f"{name}" for name, version in entries if name)


def read_run_manifest(zip_path):
    """Manifeste du builder (run_manifest.json) dans l'archive des VCF, ou {}."""
    try:
        with zipfile.ZipFile(zip_path, "r") as z:
            if RUN_MANIFEST not in z.namelist():
                return {}
            return json.loads(z.read(RUN_MANIFEST).decode("utf-8"))
    except Exception as e:
        logging.warning(f"Cannot read {RUN_MANIFEST} from '{zip_path}': {e}")
        return {}


def tgv_provenance():
    """Version et commit de TGV, empreinte de la configuration clinique."""
    return {
        "tgv_version": __version__,
        "tgv_commit": build_commit(),
        "thresholds_sha256": thresholds_sha256(),
    }
