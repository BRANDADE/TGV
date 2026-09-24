"""
Répertoire temporaire unique de la session TGV.

Tous les fichiers temporaires (exports HTML, rapport QC, IGV, graphiques) y sont
écrits ; il est supprimé à la fermeture (main.py) et, en dernier recours, via atexit.
"""
import atexit
import logging
import os
import shutil
import tempfile

_SESSION_DIR = None


def session_dir():
    """Répertoire de la session (créé à la première utilisation)."""
    global _SESSION_DIR
    if _SESSION_DIR is None or not os.path.isdir(_SESSION_DIR):
        _SESSION_DIR = tempfile.mkdtemp(prefix="tgv_")
        logging.debug(f"Session temporary directory: {_SESSION_DIR}")
    return _SESSION_DIR


def session_path(*parts):
    """Chemin sous le répertoire de session ; les répertoires parents sont créés."""
    safe_parts = [os.path.basename(str(p)) for p in parts]  # pas de '..' ni de chemin absolu
    path = os.path.join(session_dir(), *safe_parts)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def cleanup():
    """Supprime le répertoire de session et tout son contenu."""
    global _SESSION_DIR
    if _SESSION_DIR and os.path.isdir(_SESSION_DIR):
        shutil.rmtree(_SESSION_DIR, ignore_errors=True)
        logging.info(f"Session temporary directory removed: {_SESSION_DIR}")
    _SESSION_DIR = None


atexit.register(cleanup)
