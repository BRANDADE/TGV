import yaml
import logging
from functools import lru_cache
from scripts.core.config_manager import get_safe_config_path

# ============================================================
# Chargement du YAML
# ============================================================

@lru_cache(maxsize=1)  # Évite de relire et re-parser le YAML sur le disque à chaque appel sans 'data'
def load_clinical_thresholds():
    """
    Charge le fichier clinical_thresholds.yaml.
    Cherche en priorité à côté de l'exécutable (externe), puis dans les sources (interne).
    Renvoie {} si le fichier est absent, vide ou illisible.
    """

    yaml_path = get_safe_config_path("clinical_thresholds.yaml")

    try:
        with open(yaml_path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
            if data is None:
                logging.warning("clinical_thresholds.yaml is empty.")
                return {}

            logging.info(f"Clinical thresholds successfully loaded from: {yaml_path}")
            return data

    except Exception as e:
        logging.warning(f"Failed to load clinical_thresholds.yaml: {e}")
        return {}


def load_trid_aliases():
    """
    Charge configs/trid_aliases.yaml : {TRID du catalogue: bloc de clinical_thresholds.yaml}.
    Renvoie {} si le fichier est absent ou illisible (aucun alias).
    """
    path = get_safe_config_path("trid_aliases.yaml")
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}
    except FileNotFoundError:
        return {}
    except Exception as e:
        logging.warning(f"Failed to load trid_aliases.yaml: {e}")
        return {}

    if not isinstance(data, dict):
        logging.warning("trid_aliases.yaml must map catalog TRIDs to clinical_thresholds.yaml blocks.")
        return {}
    return {str(k): str(v) for k, v in data.items()}


# ============================================================
# Getter principal : bloc complet d’un TRID
# ============================================================

def get_locus_config(trid, thresholds_data=None):
    """
    Retourne le bloc YAML complet correspondant au TRID demandé.
    Exemple :
        cfg = get_locus_config("SCA1_ATXN1")

    Cette fonction :
        - charge le YAML si nécessaire
        - vérifie que le TRID existe
        - renvoie la configuration brute du locus
    """
    if thresholds_data is None:
        thresholds_data = load_clinical_thresholds()

    if trid not in thresholds_data:
        raise KeyError(f"Locus '{trid}' introuvable dans clinical_thresholds.yaml")

    return thresholds_data[trid]
