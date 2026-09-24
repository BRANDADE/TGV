"""
Version de TGV et commit exact du code exécuté (traçabilité des résultats).
"""
import os
import subprocess

__version__ = "1.1.0.dev0"

_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_COMMIT = None


def build_commit():
    """
    Commit du code exécuté :
    1. scripts/_build_info.py, généré par la CI de build (exécutable Windows) ;
    2. sinon `git rev-parse` dans le dépôt (suffixe '-dirty' si modifications locales) ;
    3. sinon 'unknown'.
    """
    global _COMMIT
    if _COMMIT is not None:
        return _COMMIT

    try:
        from scripts import _build_info  # généré par .github/workflows/build.yaml
        _COMMIT = _build_info.COMMIT
        return _COMMIT
    except Exception:
        pass

    try:
        commit = subprocess.run(
            ["git", "rev-parse", "--short=12", "HEAD"],
            cwd=_PROJECT_ROOT, capture_output=True, text=True, timeout=5, check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=_PROJECT_ROOT, capture_output=True, text=True, timeout=5, check=True,
        ).stdout.strip()
        _COMMIT = f"{commit}-dirty" if dirty else commit
    except Exception:
        _COMMIT = "unknown"
    return _COMMIT
