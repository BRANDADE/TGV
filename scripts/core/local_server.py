"""
Serveur HTTP local de TGV (IGV, graphiques TRGT), restreint à une liste blanche.

- écoute uniquement sur 127.0.0.1, port attribué par le système ;
- chaque URL commence par un jeton aléatoire propre à la session ;
- seuls les fichiers explicitement enregistrés sont servis (dictionnaire nom → chemin) :
  aucun accès au système de fichiers à partir de l'URL, aucun listing de répertoire ;
- aucun en-tête CORS : la page IGV et les BAM sont servis par la même origine ;
- requêtes Range (lecture partielle des BAM par igv.js).
"""
import http.server
import logging
import mimetypes
import os
import re
import secrets
import threading
import urllib.parse

_RANGE_RE = re.compile(r"^bytes=(\d*)-(\d*)$")
_CHUNK = 64 * 1024


class _Handler(http.server.BaseHTTPRequestHandler):
    server_version = "TGV"
    sys_version = ""

    def do_GET(self):
        self._serve(send_body=True)

    def do_HEAD(self):
        self._serve(send_body=False)

    def log_message(self, fmt, *args):
        logging.debug("Local server: " + fmt % args)

    def _resolve(self):
        path = urllib.parse.urlsplit(self.path).path
        prefix = f"/{self.server.token}/"
        if not path.startswith(prefix):
            return None
        name = urllib.parse.unquote(path[len(prefix):])
        file_path = self.server.files.get(name)
        if file_path is None or not os.path.isfile(file_path):
            return None
        return file_path

    def _serve(self, send_body):
        file_path = self._resolve()
        if file_path is None:
            self.send_error(404, "Not found")
            return

        size = os.path.getsize(file_path)
        start, end = 0, size - 1
        status = 200

        range_header = self.headers.get("Range")
        if range_header:
            match = _RANGE_RE.match(range_header.strip())
            if not match or (not match.group(1) and not match.group(2)):
                self.send_error(416, "Invalid range")
                return
            first, last = match.groups()
            if first:
                start = int(first)
                end = min(int(last), size - 1) if last else size - 1
            else:  # suffixe : les N derniers octets
                start = max(0, size - int(last))
            if start >= size or start > end:
                self.send_response(416)
                self.send_header("Content-Range", f"bytes */{size}")
                self.end_headers()
                return
            status = 206

        ctype = mimetypes.guess_type(file_path)[0] or "application/octet-stream"
        self.send_response(status)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(end - start + 1 if size else 0))
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Cache-Control", "no-store")
        if status == 206:
            self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
        self.end_headers()

        if not send_body or size == 0:
            return
        with open(file_path, "rb") as f:
            f.seek(start)
            remaining = end - start + 1
            while remaining > 0:
                chunk = f.read(min(_CHUNK, remaining))
                if not chunk:
                    break
                self.wfile.write(chunk)
                remaining -= len(chunk)


class LocalFileServer:
    """Serveur de fichiers en liste blanche (voir docstring du module)."""

    def __init__(self):
        self._httpd = None
        self._lock = threading.Lock()
        self.token = secrets.token_urlsafe(24)
        self.files = {}

    @property
    def port(self):
        return self._httpd.server_address[1] if self._httpd else None

    def start(self):
        with self._lock:
            if self._httpd is None:
                httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), _Handler)
                httpd.daemon_threads = True
                httpd.token = self.token
                httpd.files = self.files
                threading.Thread(target=httpd.serve_forever, kwargs={"poll_interval": 0.1}, daemon=True).start()
                self._httpd = httpd
                logging.info(f"Local file server started on 127.0.0.1:{self.port}")
        return self

    def register(self, name, file_path):
        """Autorise le fichier sous le nom donné et renvoie son URL."""
        self.start()
        self.files[name] = os.path.abspath(file_path)
        return self.url(name)

    def unregister_prefix(self, prefix):
        for name in [n for n in self.files if n.startswith(prefix)]:
            del self.files[name]

    def url(self, name):
        return f"http://127.0.0.1:{self.port}/{self.token}/{urllib.parse.quote(name)}"

    def stop(self):
        with self._lock:
            if self._httpd is not None:
                try:
                    self._httpd.shutdown()
                    self._httpd.server_close()
                except Exception:
                    pass
                self._httpd = None
            self.files.clear()


_SERVER = None


def get_server():
    """Serveur unique de la session, démarré à la première utilisation."""
    global _SERVER
    if _SERVER is None:
        _SERVER = LocalFileServer()
    return _SERVER.start()


def stop_server():
    global _SERVER
    if _SERVER is not None:
        _SERVER.stop()
        _SERVER = None
