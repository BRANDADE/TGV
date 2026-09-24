import os
import io
import time
import zipfile
import logging
import webbrowser

from scripts.core.local_server import get_server
from scripts.core.session_tmp import session_path

# Préfixe des graphiques dans la liste blanche du serveur local
PLOTS_PREFIX = "plots/"


def plots_prefix(sample_name):
    return f"{PLOTS_PREFIX}{os.path.basename(sample_name)}/"


def open_svg(zip_path, inner_zip, svg_file, sample_name):
    """
    Extrait un graphique TRGT dans le répertoire de session et l'ouvre via le
    serveur local (liste blanche : seul ce fichier est ajouté).
    """
    logging.info(f"Extracting SVG plot '{svg_file}' from nested archive: '{inner_zip}'")
    try:
        with zipfile.ZipFile(zip_path, "r") as outer:
            with zipfile.ZipFile(io.BytesIO(outer.read(inner_zip))) as inner:
                svg_bytes = inner.read(svg_file)

        svg_tmp_path = session_path("plots", sample_name, svg_file)
        logging.debug(f"Writing temporary SVG plot file to: {svg_tmp_path}")
        with open(svg_tmp_path, "wb") as f:
            f.write(svg_bytes)

        url = get_server().register(plots_prefix(sample_name) + os.path.basename(svg_file), svg_tmp_path)
        logging.info(f"Launching web browser for SVG visualization at: {url}")
        webbrowser.open(f"{url}?t={time.time()}")

    except Exception as e:
        logging.error(f"Failed to extract or serve SVG plot '{svg_file}': {e}", exc_info=True)


def forget_plots(sample_name):
    """Retire les graphiques d'un patient du serveur et du répertoire de session."""
    get_server().unregister_prefix(plots_prefix(sample_name))
    directory = os.path.dirname(session_path("plots", sample_name, "x"))
    for name in os.listdir(directory):
        try:
            os.remove(os.path.join(directory, name))
        except OSError:
            pass
