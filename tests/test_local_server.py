"""D1/D2 — serveur local en liste blanche et répertoire temporaire de session."""
import http.client
import os

import pytest

from scripts.core import session_tmp
from scripts.core.local_server import LocalFileServer


@pytest.fixture
def server(tmp_path):
    srv = LocalFileServer().start()
    data = tmp_path / "S1.trgt.spanning.sorted.bam"
    data.write_bytes(bytes(range(256)) * 4)
    srv.register("igv/S1.trgt.spanning.sorted.bam", str(data))
    (tmp_path / "secret.txt").write_text("secret")
    yield srv
    srv.stop()


def _get(srv, path, headers=None):
    conn = http.client.HTTPConnection("127.0.0.1", srv.port, timeout=5)
    conn.request("GET", path, headers=headers or {})
    resp = conn.getresponse()
    body = resp.read()
    conn.close()
    return resp, body


def test_registered_file_is_served_without_cors(server):
    resp, body = _get(server, f"/{server.token}/igv/S1.trgt.spanning.sorted.bam")
    assert resp.status == 200
    assert len(body) == 1024
    assert resp.getheader("Access-Control-Allow-Origin") is None


def test_range_request(server):
    resp, body = _get(server, f"/{server.token}/igv/S1.trgt.spanning.sorted.bam", {"Range": "bytes=10-19"})
    assert resp.status == 206
    assert body == bytes(range(10, 20))
    assert resp.getheader("Content-Range") == "bytes 10-19/1024"


def test_suffix_range_request(server):
    resp, body = _get(server, f"/{server.token}/igv/S1.trgt.spanning.sorted.bam", {"Range": "bytes=-4"})
    assert resp.status == 206
    assert body == bytes([252, 253, 254, 255])


@pytest.mark.parametrize("path", [
    "/{token}/..%2f..%2f..%2f..%2fetc/os-release",
    "/{token}/../secret.txt",
    "/{token}/igv/",
    "/{token}/",
    "/genome/..%2f..%2f..%2f..%2fetc/os-release",
    "/igv/S1.trgt.spanning.sorted.bam",          # sans jeton
    "/wrongtoken/igv/S1.trgt.spanning.sorted.bam",
])
def test_everything_else_is_404(server, path):
    resp, _ = _get(server, path.format(token=server.token))
    assert resp.status == 404


def test_unregister_prefix(server):
    server.unregister_prefix("igv/")
    resp, _ = _get(server, f"/{server.token}/igv/S1.trgt.spanning.sorted.bam")
    assert resp.status == 404


def test_session_directory_is_removed_on_cleanup():
    path = session_tmp.session_path("exports", "tgv_export_S1.html")
    with open(path, "w") as f:
        f.write("x")
    root = session_tmp.session_dir()
    assert os.path.basename(root).startswith("tgv_")
    session_tmp.cleanup()
    assert not os.path.exists(root)


def test_session_path_cannot_escape():
    path = session_tmp.session_path("..", "..", "etc", "passwd")
    assert os.path.commonpath([path, session_tmp.session_dir()]) == session_tmp.session_dir()
    session_tmp.cleanup()
