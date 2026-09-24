# main.py
import sys
import os
import logging
import scripts.core.i18n  # noqa: F401  (traductions et bouton de langue appliqués à PySimpleGUI)

# Importation du configurateur de logs
from scripts.core.logger import setup_logging
from scripts.ui.main_window import run_main_window


# ==============================================================================
# OPTION DE DEBUGGING WINDOWS : ALLOCATION DE CONSOLE
# ==============================================================================
# Ce bloc de code permet d'ouvrir de force une invite de commande (console) sous 
# Windows, même si l'application est compilée en mode "Sans Console" (ex: PyInstaller -w).
#
# USAGE POUR LE DIAGNOSTIC DES BUGS EN PHASE DE TEST :
# 1. Décommentez les lignes d'imports, la fonction 'open_console' et son appel.
# 2. Lancez l'application. Une console Windows s'ouvrira en parallèle de l'IHM
#    et affichera tous les flux de logs, de print() et d'erreurs système en temps réel.
# ==============================================================================
# import ctypes
# import sys

# def open_console():
#     """
#     Alloue dynamiquement une console système sous Windows pour afficher 
#     en temps réel les logs et les erreurs de l'application (résolution de bugs).
#     """
#     # Appel à l'API Windows (kernel32.dll) pour créer et allouer une console active
#     ctypes.windll.kernel32.AllocConsole()
#     # Redirection des flux d'écriture standards de Python vers la console créée (CONOUT$)
#     sys.stdout = open("CONOUT$", "w", encoding="utf-8")
#     sys.stderr = open("CONOUT$", "w", encoding="utf-8")

# # Décommentez l'appel ci-dessous pour activer la console de diagnostic sous Windows
# # open_console()
# ==============================================================================


# --- REDIRECTION DES LOGS POUR PYINSTALLER (SANS CONSOLE) ---
# Empêche l'application de crasher sur Windows lors de l'écriture si aucun terminal n'est actif
if sys.stdout is None:
    sys.stdout = open(os.devnull, "w")
if sys.stderr is None:
    sys.stderr = open(os.devnull, "w")


def cleanup_temp_files():
    """
    Arrête le serveur local et supprime le répertoire temporaire de la session
    (exports HTML, rapport QC, IGV, graphiques). Les échecs sont consignés sans
    interrompre la fermeture.
    """
    logging.info("Starting cleanup of temporary files...")
    try:
        from scripts.core.local_server import stop_server
        stop_server()
    except Exception as e:
        logging.warning(f"Failed to stop local server: {e}")
    try:
        from scripts.core import session_tmp
        session_tmp.cleanup()
    except Exception as e:
        logging.warning(f"Failed to remove session temporary directory: {e}")


def main():
    # Initialisation globale du système de logging quotidien
    setup_logging()
    
    try:
        logging.info("Launching main graphical user interface...")
        run_main_window()
        logging.info("Application closed normally.")
    except Exception:
        # Enregistrement du crash complet dans l'unique fichier de log de l'application
        logging.critical("A critical error occurred during execution:", exc_info=True)
        sys.exit(1)
    finally:
        # Exécuté dans tous les cas à la fermeture de l'outil
        cleanup_temp_files()


if __name__ == "__main__":
    main()
