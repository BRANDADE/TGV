import os
import sys
import zipfile
import tempfile
import shutil
import webbrowser
import threading
import http.server
import socketserver
import re
import json
import socket
import logging
import urllib.parse
import PySimpleGUI as sg
from scripts.core.i18n import tr

PADDING = 75
CURRENT_TMPDIR = None
CURRENT_SERVER = None
CURRENT_GENOME_DIR = None


def is_online():
    logging.debug("Probing active internet connection...")
    try:
        socket.setdefaulttimeout(1.5)
        socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect(("s3.amazonaws.com", 443))
        return True
    except socket.error:
        return False


def get_val(obj, attr_name, default="N/A"):
    if not obj:
        return default
    val = getattr(obj, attr_name, None)
    if val is None and isinstance(obj, dict):
        val = obj.get(attr_name)
    if val in (None, "", ".", "None", "N/A"):
        return default
    return str(val)


class RangeRequestHandler(http.server.SimpleHTTPRequestHandler):
    def send_head(self):
        path = self.translate_path(self.path)
        if not os.path.exists(path) or os.path.isdir(path):
            return super().send_head()
            
        ctype = self.guess_type(path)
        try:
            f = open(path, 'rb')
        except OSError:
            self.send_error(404, "File not found")
            return None
            
        fs = os.fstat(f.fileno())
        file_size = fs[6]
        
        range_header = self.headers.get('Range')
        if range_header:
            match = re.match(r'bytes=(\d*)-(\d*)', range_header)
            if match:
                start, end = match.groups()
                start = int(start) if start else 0
                end = int(end) if end else file_size - 1
                end = min(end, file_size - 1)
                
                if start >= file_size or start > end:
                    self.send_error(416, "Requested range not satisfiable")
                    f.close()
                    return None
                    
                self.send_response(206)
                self.send_header('Content-Type', ctype)
                self.send_header('Accept-Ranges', 'bytes')
                self.send_header('Content-Range', f'bytes {start}-{end}/{file_size}')
                self.send_header('Content-Length', str(end - start + 1))
                self.send_header('Access-Control-Allow-Origin', '*')
                self.end_headers()
                
                f.seek(start)
                self.range_info = (start, end)
                return f
                
        self.send_response(200)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(file_size))
        self.send_header('Accept-Ranges', 'bytes')
        self.send_header('Access-Control-Allow-Origin', '*')
        self.end_headers()
        self.range_info = None
        return f

    def copyfile(self, source, outputfile):
        if hasattr(self, 'range_info') and self.range_info:
            start, end = self.range_info
            remaining = end - start + 1
            buffer_size = 64 * 1024
            while remaining > 0:
                chunk = source.read(min(remaining, buffer_size))
                if not chunk:
                    break
                outputfile.write(chunk)
                remaining -= len(chunk)
        else:
            super().copyfile(source, outputfile)

    def translate_path(self, path):
        global CURRENT_GENOME_DIR
        decoded_path = urllib.parse.unquote(path)
        if decoded_path.startswith("/genome/") and CURRENT_GENOME_DIR:
            relative_path = decoded_path[len("/genome/"):]
            return os.path.join(CURRENT_GENOME_DIR, relative_path)
        return super().translate_path(path)


def start_local_server(directory):
    global CURRENT_SERVER
    stop_local_server()
    handler = lambda *args, **kwargs: RangeRequestHandler(*args, directory=directory, **kwargs)
    try:
        server = socketserver.TCPServer(("127.0.0.1", 0), handler)
        port = server.socket.getsockname()[1]
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        CURRENT_SERVER = server
        return port
    except Exception as e:
        logging.error(f"Failed to start local HTTP server: {e}")
        raise e


def stop_local_server():
    global CURRENT_SERVER
    if CURRENT_SERVER:
        try:
            CURRENT_SERVER.shutdown()
            CURRENT_SERVER.server_close()
        except Exception:
            pass
        CURRENT_SERVER = None


def cleanup_tmpdir_force():
    global CURRENT_TMPDIR
    stop_local_server()
    if CURRENT_TMPDIR:
        shutil.rmtree(CURRENT_TMPDIR, ignore_errors=True)
        CURRENT_TMPDIR = None


def get_asset_path(filename):
    if hasattr(sys, '_MEIPASS'):
        return os.path.join(sys._MEIPASS, 'assets', filename)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    local_assets_dir = os.path.join(script_dir, 'assets', filename)
    if os.path.exists(local_assets_dir):
        return local_assets_dir
    local_direct = os.path.join(script_dir, filename)
    if os.path.exists(local_direct):
        return local_direct
    dev_root_dir = os.path.abspath(os.path.join(script_dir, "..", ".."))
    return os.path.join(dev_root_dir, 'assets', filename)


def get_igv_js_version():
    asset_path = get_asset_path("igv.min.js")
    if os.path.exists(asset_path):
        try:
            with open(asset_path, "r", encoding="utf-8", errors="ignore") as f:
                header = f.read(500)
                match = re.search(r"igv\.js\s+v?([0-9.]+)", header)
                if match:
                    return match.group(1)
        except Exception:
            pass
    return "unknown"


def open_igv(genome_fasta_path=None, 
             spanning_zip_path=None, spanning_bam_file=None, spanning_bai_file=None, 
             mapped_zip_path=None, mapped_bam_file=None, mapped_bai_file=None, 
             chrom=None, start=None, end=None,
             sample_name=None,
             row=None):
             
    has_local_fasta = False
    fasta_filename = None
    fai_filename = None
    
    if genome_fasta_path and os.path.exists(genome_fasta_path):
        fasta_filename = os.path.basename(genome_fasta_path)
        fai_filename = fasta_filename + ".fai"
        fai_full_path = os.path.join(os.path.dirname(genome_fasta_path), fai_filename)
        if os.path.exists(fai_full_path):
            has_local_fasta = True

    online_mode = is_online()

    if not has_local_fasta and not online_mode:
        sg.popup(tr("Erreur"), tr("Veuillez sélectionner un fichier génome de référence local."))
        return

    r_obj = None
    if row and isinstance(row, dict):
        r_obj = row.get("Result_obj")
        if r_obj:
            if not chrom: chrom = getattr(r_obj, "chrom", None)
            if not start: start = getattr(r_obj, "start", None)
            if not end: end = getattr(r_obj, "end", None)
            if not sample_name: 
                sample_name = getattr(r_obj, "sample_id", None) or getattr(r_obj, "sample_name", None)

    if chrom is None or start is None or end is None:
        sg.popup(tr("Aucun locus sélectionné"))
        return
        
    global CURRENT_TMPDIR, CURRENT_GENOME_DIR

    if not sample_name:
        if spanning_bam_file:
            sample_name = os.path.basename(spanning_bam_file).split('.')[0]
        elif mapped_bam_file:
            sample_name = os.path.basename(mapped_bam_file).split('.')[0]
        else:
            sample_name = "Sample"

    if has_local_fasta:
        CURRENT_GENOME_DIR = os.path.dirname(genome_fasta_path)
        genome_option_js = json.dumps({
            "id": "local_genome",
            "name": fasta_filename,
            "fastaURL": f"/genome/{fasta_filename}",
            "indexURL": f"/genome/{fai_filename}"
        }, indent=4)
    else:
        genome_option_js = '"hg38"'

    cleanup_tmpdir_force()
    CURRENT_TMPDIR = tempfile.mkdtemp()

    asset_js_path = get_asset_path("igv.min.js")
    if not os.path.exists(asset_js_path):
        sg.popup(tr("Erreur"), f"igv.min.js introuvable : {asset_js_path}")
        return

    try:
        shutil.copy(asset_js_path, os.path.join(CURRENT_TMPDIR, "igv.min.js"))
    except Exception as e:
        sg.popup(f"Erreur :\n{e}")
        return

    tracks = []
    if online_mode:
        tracks.append({
            "name": "Reference Genes (UCSC RefSeq)",
            "type": "annotation",
            "format": "refgene",
            "url": "https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefSeq.sorted.txt.gz",
            "indexURL": "https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefSeq.sorted.txt.gz.tbi",
            "displayMode": "COLLAPSED",
            "order": 1
        })

    if spanning_zip_path and spanning_bam_file and spanning_bai_file:
        try:
            with zipfile.ZipFile(spanning_zip_path, "r") as z:
                z.extract(spanning_bam_file, CURRENT_TMPDIR)
                z.extract(spanning_bai_file, CURRENT_TMPDIR)
            tracks.append({
                "name": f"Spanning BAM - {os.path.basename(spanning_bam_file)}",
                "url": f"./{spanning_bam_file}",
                "indexURL": f"./{spanning_bai_file}",
                "type": "alignment",
                "format": "bam",
                "colorBy": "strand",
                "displayMode": "EXPANDED",
                "alignmentRowHeight": 14,
                "autoHeight": True,
                "maxHeight": 600,
                "samplingDepth": 1000,
                "showSoftClips": True,
                "showDeletionText": True,
                "showInsertionText": True,
                "sort": {"option": "INSERT_SIZE", "direction": "ASC"}
            })
        except Exception as e:
            sg.popup(f"Erreur extraction spanning BAM :\n{e}")

    if mapped_zip_path and mapped_bam_file and mapped_bai_file:
        try:
            with zipfile.ZipFile(mapped_zip_path, "r") as z:
                z.extract(mapped_bam_file, CURRENT_TMPDIR)
                z.extract(mapped_bai_file, CURRENT_TMPDIR)
            tracks.append({
                "name": f"Mapped BAM - {os.path.basename(mapped_bam_file)}",
                "url": f"./{mapped_bam_file}",
                "indexURL": f"./{mapped_bai_file}",
                "type": "alignment",
                "format": "bam",
                "colorBy": "strand",
                "displayMode": "EXPANDED",
                "alignmentRowHeight": 14,
                "autoHeight": True,
                "maxHeight": 600,
                "samplingDepth": 1000,
                "showSoftClips": True,
                "showDeletionText": True,
                "showInsertionText": True,
                "sort": {"option": "INSERT_SIZE", "direction": "ASC"}
            })
        except Exception as e:
            sg.popup(f"Erreur extraction mapped BAM :\n{e}")

    if not tracks:
        sg.popup(tr("Aucun BAM disponible"))
        return

    # Noms réels des fichiers affichés : l'identité du patient doit être vérifiable
    shown_files = [os.path.basename(f) for f in (spanning_bam_file, mapped_bam_file) if f]
    files_str = ", ".join(shown_files) if shown_files else "N/A"

    start_padded = max(0, start - PADDING)
    end_padded = end + PADDING
    region = f"{chrom}:{start_padded}-{end_padded}"

    # --- Construction du Dashboard d'en-tête traduit ---
    if row and isinstance(row, dict):
        locus_title = row.get("Locus") or f"{chrom}:{start}-{end}"
        motifs_str = row.get("Motifs") or "N/A"
        
        rep1 = row.get("Allèle 1 - Répétition") or "N/A"
        if rep1 in (None, "", "None", "N/A") and r_obj:
            rep1 = get_val(r_obj, "rep1_raw")
            
        depth1 = "N/A"
        meth1_pct = "N/A"
        purity1_pct = "N/A"
        
        applied_class = row.get("Classification", "")
        if applied_class and isinstance(applied_class, str) and "/" in applied_class:
            class1, class2 = [x.strip() for x in applied_class.split("/")]
        else:
            class1 = str(applied_class) if applied_class else tr("Non classifié")
            class2 = tr("Non classifié")
            
        applied_gt = row.get("Génotype", "")
        if applied_gt and isinstance(applied_gt, str) and "/" in applied_gt:
            gen1, gen2 = [x.strip() for x in applied_gt.split("/")]
        else:
            gen1 = str(applied_gt) if applied_gt else "None"
            gen2 = "None"

        rep2 = row.get("Allèle 2 - Répétition") or "N/A"
        if rep2 in (None, "", "None", "N/A") and r_obj:
            rep2 = get_val(r_obj, "rep2_raw")
            
        depth2 = "N/A"
        meth2_pct = "N/A"
        purity2_pct = "N/A"

        if r_obj:
            depth1 = get_val(r_obj, "depth1_raw")
            depth2 = get_val(r_obj, "depth2_raw")
            
            meth1 = get_val(r_obj, "methylation1")
            if meth1 != "N/A":
                try: meth1_pct = f"{int(float(meth1) * 100)}%"
                except ValueError: meth1_pct = meth1
                
            meth2 = get_val(r_obj, "methylation2")
            if meth2 != "N/A":
                try: meth2_pct = f"{int(float(meth2) * 100)}%"
                except ValueError: meth2_pct = meth2
                
            pur1 = get_val(r_obj, "purity1")
            if pur1 != "N/A":
                try: purity1_pct = f"{int(float(pur1) * 100)}%"
                except ValueError: purity1_pct = pur1
                
            pur2 = get_val(r_obj, "purity2")
            if pur2 != "N/A":
                try: purity2_pct = f"{int(float(pur2) * 100)}%"
                except ValueError: purity2_pct = pur2
            
            if class1 in (None, "", "None", "Non classé", "Non classifié"):
                class1 = get_val(r_obj, "classification1_bio") or get_val(r_obj, "classification1_raw") or tr("Non classifié")
            if class2 in (None, "", "None", "Non classé", "Non classifié"):
                class2 = get_val(r_obj, "classification2_bio") or get_val(r_obj, "classification2_raw") or tr("Non classifié")
                
            if gen1 in (None, "", "None"):
                gen1 = get_val(r_obj, "genotype1_bio") or get_val(r_obj, "genotype1_raw") or "None"
            if gen2 in (None, "", "None"):
                gen2 = get_val(r_obj, "genotype2_bio") or get_val(r_obj, "genotype2_raw") or "None"

        if class1 in (None, "", "None", "Non classifié", "Non classé"):
            class1_html = f'<span style="font-size: 0.72rem; background: #475569; color: #cbd5e1; padding: 2px 8px; border-radius: 9999px; font-weight: 700;">{tr("Non classifié")}</span>'
        else:
            class1_color = "#f43f5e" if "patho" in class1.lower() else "#10b981" if any(w in class1.lower() for w in ["benign", "sain"]) else "#f59e0b"
            class1_html = f'<span style="font-size: 0.72rem; background: {class1_color}; color: white; padding: 2px 8px; border-radius: 9999px; font-weight: 700;">{tr(class1)}</span>'
            
        gen1_html = ""
        if gen1 not in (None, "", "None"):
            gen1_html = f'<div style="font-size: 0.8rem; color: #cbd5e1; margin-top: 4px; font-weight: 500;">{tr("Génotype :")} <strong style="color: #38bdf8;">{gen1}</strong></div>'

        if class2 in (None, "", "None", "Non classifié", "Non classé"):
            class2_html = f'<span style="font-size: 0.72rem; background: #475569; color: #cbd5e1; padding: 2px 8px; border-radius: 9999px; font-weight: 700;">{tr("Non classifié")}</span>'
        else:
            class2_color = "#f43f5e" if "patho" in class2.lower() else "#10b981" if any(w in class2.lower() for w in ["benign", "sain"]) else "#f59e0b"
            class2_html = f'<span style="font-size: 0.72rem; background: {class2_color}; color: white; padding: 2px 8px; border-radius: 9999px; font-weight: 700;">{tr(class2)}</span>'
            
        gen2_html = ""
        if gen2 not in (None, "", "None"):
            gen2_html = f'<div style="font-size: 0.8rem; color: #cbd5e1; margin-top: 4px; font-weight: 500;">{tr("Génotype :")} <strong style="color: #38bdf8;">{gen2}</strong></div>'

        metrics1_html = f'<div>{tr("Profondeur :")} <strong style="color: #f8fafc;">{depth1}x</strong></div>' if depth1 != "N/A" else f'<div>{tr("Profondeur :")} <strong style="color: #f8fafc;">N/A</strong></div>'
        if meth1_pct != "N/A":
            metrics1_html += f'<div style="margin-left: 12px; border-left: 1px solid rgba(255,255,255,0.15); padding-left: 12px;">{tr("Méthylation :")} <strong style="color: #34d399;">{meth1_pct}</strong></div>'
        if purity1_pct != "N/A":
            metrics1_html += f'<div style="margin-left: 12px; border-left: 1px solid rgba(255,255,255,0.15); padding-left: 12px;">{tr("Pureté :")} <strong style="color: #67e8f9;">{purity1_pct}</strong></div>'

        metrics2_html = f'<div>{tr("Profondeur :")} <strong style="color: #f8fafc;">{depth2}x</strong></div>' if depth2 != "N/A" else f'<div>{tr("Profondeur :")} <strong style="color: #f8fafc;">N/A</strong></div>'
        if meth2_pct != "N/A":
            metrics2_html += f'<div style="margin-left: 12px; border-left: 1px solid rgba(255,255,255,0.15); padding-left: 12px;">{tr("Méthylation :")} <strong style="color: #34d399;">{meth2_pct}</strong></div>'
        if purity2_pct != "N/A":
            metrics2_html += f'<div style="margin-left: 12px; border-left: 1px solid rgba(255,255,255,0.15); padding-left: 12px;">{tr("Pureté :")} <strong style="color: #67e8f9;">{purity2_pct}</strong></div>'

        clinical_header_html = f"""
        <div style="background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%); color: white; padding: 20px 25px; border-radius: 8px; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1); margin-bottom: 20px; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;">
            <div style="display: flex; justify-content: space-between; align-items: center; border-bottom: 1px solid rgba(255, 255, 255, 0.15); padding-bottom: 15px; margin-bottom: 15px;">
                <div style="display: flex; align-items: center; gap: 16px;">
                    <div style="width: 36px; height: 36px; border-radius: 50%; background: linear-gradient(135deg, #e8457a 0%, #9e1f4f 100%); flex-shrink: 0;"></div>
                    <div>
                        <h1 style="margin: 0; font-size: 1.5rem; font-weight: 700; color: #f8fafc; letter-spacing: -0.025em;">TGV - TRGT Global Viewer</h1>
                        <p style="margin: 4px 0 0 0; font-size: 0.85rem; color: #94a3b8;">{tr("Visualisation IGV")}</p>
                    </div>
                </div>
                <div style="text-align: right;">
                    <div style="font-size: 1.25rem; font-weight: 700; color: #e8457a;">{tr("Patient :")} <span style="color: #e8457a;">{sample_name}</span></div>
                    <div style="font-size: 0.85rem; color: #cbd5e1; margin-top: 4px; font-weight: 500;">
                        {tr("Locus :")} <span style="font-family: monospace; background: rgba(255,255,255,0.15); padding: 2px 6px; border-radius: 4px;">{chrom}:{start}-{end}</span>
                    </div>
                    <div style="font-size: 0.8rem; color: #cbd5e1; margin-top: 4px;">
                        {tr("Fichiers :")} <span style="font-family: monospace;">{files_str}</span>
                    </div>
                </div>
            </div>
            
            <div style="display: grid; grid-template-columns: 1fr 1.5fr 1.5fr; gap: 20px;">
                <div style="background: rgba(255, 255, 255, 0.03); border: 1px solid rgba(255, 255, 255, 0.1); border-radius: 6px; padding: 12px 15px;">
                    <div style="font-size: 0.72rem; text-transform: uppercase; color: #94a3b8; font-weight: 700; letter-spacing: 0.05em; margin-bottom: 8px;">{tr("MÉTADONNÉES LOCUS")}</div>
                    <div style="font-size: 1.15rem; font-weight: 700; color: #f8fafc; margin-bottom: 6px;">{locus_title}</div>
                    <div style="font-size: 0.85rem; color: #cbd5e1; margin-top: 4px;">{tr("Motif :")} <span style="font-family: monospace; background: rgba(255,255,255,0.1); padding: 1px 4px; border-radius: 3px; color: #67e8f9;">{motifs_str}</span></div>
                </div>
                
                <div style="background: rgba(255, 255, 255, 0.03); border: 1px solid rgba(255, 255, 255, 0.1); border-radius: 6px; padding: 12px 15px; display: flex; flex-direction: column; justify-content: space-between;">
                    <div>
                        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;">
                            <span style="font-size: 0.72rem; text-transform: uppercase; color: #94a3b8; font-weight: 700; letter-spacing: 0.05em;">{tr("ALLÈLE 1")}</span>
                            {class1_html}
                        </div>
                        <div style="font-size: 1.35rem; font-weight: 800; color: #38bdf8; font-family: monospace;">{rep1}</div>
                        {gen1_html}
                    </div>
                    <div style="display: flex; font-size: 0.8rem; color: #cbd5e1; margin-top: 8px; border-top: 1px solid rgba(255,255,255,0.08); padding-top: 6px;">
                        {metrics1_html}
                    </div>
                </div>

                <div style="background: rgba(255, 255, 255, 0.03); border: 1px solid rgba(255, 255, 255, 0.1); border-radius: 6px; padding: 12px 15px; display: flex; flex-direction: column; justify-content: space-between;">
                    <div>
                        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;">
                            <span style="font-size: 0.72rem; text-transform: uppercase; color: #94a3b8; font-weight: 700; letter-spacing: 0.05em;">{tr("ALLÈLE 2")}</span>
                            {class2_html}
                        </div>
                        <div style="font-size: 1.35rem; font-weight: 800; color: #38bdf8; font-family: monospace;">{rep2}</div>
                        {gen2_html}
                    </div>
                    <div style="display: flex; font-size: 0.8rem; color: #cbd5e1; margin-top: 8px; border-top: 1px solid rgba(255,255,255,0.08); padding-top: 6px;">
                        {metrics2_html}
                    </div>
                </div>
            </div>
        </div>
        """
    else:
        clinical_header_html = f"""
        <div style="background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%); color: white; padding: 18px 25px; border-radius: 8px; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1); margin-bottom: 20px; display: flex; justify-content: space-between; align-items: center; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;">
            <div style="display: flex; align-items: center; gap: 16px;">
                <div style="width: 36px; height: 36px; border-radius: 50%; background: linear-gradient(135deg, #e8457a 0%, #9e1f4f 100%); flex-shrink: 0;"></div>
                <div>
                    <h1 style="margin: 0; font-size: 1.5rem; font-weight: 700; color: #f8fafc;">TGV - TRGT Global Viewer</h1>
                    <p style="margin: 4px 0 0 0; font-size: 0.85rem; color: #94a3b8;">{tr("Visualisation IGV")}</p>
                </div>
            </div>
            <div style="text-align: right;">
                <div style="font-size: 1.2rem; font-weight: 700; color: #e8457a;">{tr("ID Patient :")} <span style="color: #e8457a;">{sample_name}</span></div>
                <div style="font-size: 0.85rem; color: #cbd5e1; margin-top: 4px; font-weight: 500;">
                    {tr("Locus :")} <span style="font-family: monospace; background: rgba(255,255,255,0.15); padding: 2px 6px; border-radius: 4px;">{chrom}:{start}-{end}</span>
                </div>
                <div style="font-size: 0.8rem; color: #cbd5e1; margin-top: 4px;">
                    {tr("Fichiers :")} <span style="font-family: monospace;">{files_str}</span>
                </div>
            </div>
        </div>
        """

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>{tr("Visualisation IGV")} - {sample_name}</title>
    <script src="./igv.min.js"></script>
</head>
<body style="margin:0; padding:15px; background-color: #f4f6f9; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;">
    {clinical_header_html}
    <div id="igv-div" style="background-color: white; padding: 15px; border-radius: 8px; box-shadow: 0 10px 15px -3px rgba(0,0,0,0.05); border: 1px solid #e2e8f0;"></div>

    <script type="text/javascript">
        document.addEventListener("DOMContentLoaded", function () {{
            var div = document.getElementById("igv-div");
            var options = {{
                genome: {genome_option_js},
                locus: "{region}",
                tracks: {json.dumps(tracks, indent=4)}
            }};

            igv.createBrowser(div, options);
        }});
    </script>
</body>
</html>
"""

    with open(os.path.join(CURRENT_TMPDIR, "index.html"), "w", encoding="utf-8") as f:
        f.write(html_content)

    try:
        port = start_local_server(CURRENT_TMPDIR)
        url = f"http://127.0.0.1:{port}/index.html"
        webbrowser.open(url)
    except Exception as e:
        sg.popup(f"Impossible de lancer igv.js :\n{e}")