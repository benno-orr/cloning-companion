"""Local, transactional storage for completed design runs (not exported files)."""
from __future__ import annotations

import base64
import copy
import io
import json
import sqlite3
import uuid
import zipfile
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path


@contextmanager
def _database(root: Path):
    folder = root / "DesignLibrary"
    folder.mkdir(parents=True, exist_ok=True, mode=0o700)
    path = folder / "designs.sqlite3"
    connection = sqlite3.connect(path)
    try:
        path.chmod(0o600)
        with connection:
            connection.execute("""CREATE TABLE IF NOT EXISTS runs (
                id TEXT PRIMARY KEY, created TEXT NOT NULL, title TEXT NOT NULL,
                plasmid_count INTEGER NOT NULL, result_json TEXT NOT NULL, files_zip BLOB NOT NULL
            )""")
            yield connection
    finally:
        connection.close()


def save(root: Path, result: dict, files: dict) -> dict:
    run = {"id": uuid.uuid4().hex, "created": datetime.now(timezone.utc).isoformat(),
           "title": result.get("project", "Golden Gate design"), "plasmidCount": result.get("plasmidCount", 0)}
    stored = copy.deepcopy(result)
    stored.update(savedRun=run, temporary=False)
    for graphic in stored.get("graphics", []):
        graphic.pop("src", None)  # Reconstructed from the archived PNG bytes.
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for key, item in files.items():
            archive.writestr(key, item["data"])
    with _database(root) as db:
        db.execute("INSERT INTO runs VALUES (?, ?, ?, ?, ?, ?)",
                   (run["id"], run["created"], run["title"], run["plasmidCount"], json.dumps(stored), buffer.getvalue()))
    return run


def list_runs(root: Path) -> list[dict]:
    if not (root / "DesignLibrary" / "designs.sqlite3").is_file():
        return []
    with _database(root) as db:
        return [{"id": row[0], "created": row[1], "title": row[2], "plasmidCount": row[3]}
                for row in db.execute("SELECT id, created, title, plasmid_count FROM runs ORDER BY created DESC, id DESC")]


def load(root: Path, run_id: str) -> tuple[dict, dict]:
    if not (root / "DesignLibrary" / "designs.sqlite3").is_file():
        raise ValueError("No saved designs yet")
    with _database(root) as db:
        row = db.execute("SELECT result_json, files_zip FROM runs WHERE id = ?", (run_id,)).fetchone()
    if row is None:
        raise ValueError("Saved design not found")
    result = json.loads(row[0])
    with zipfile.ZipFile(io.BytesIO(row[1])) as archive:
        files = {info["id"]: {**info, "data": archive.read(info["id"])} for info in result["downloads"]}
    for graphic in result.get("graphics", []):
        item = files[graphic["downloadId"]]
        mime = "image/svg+xml" if item["filename"].endswith(".svg") else "image/png"
        graphic["src"] = f"data:{mime};base64," + base64.b64encode(item["data"]).decode("ascii")
    return result, files
