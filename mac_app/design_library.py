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
            if "project_id" not in {row[1] for row in connection.execute("PRAGMA table_info(runs)")}:
                connection.execute("ALTER TABLE runs ADD COLUMN project_id TEXT")
            yield connection
    finally:
        connection.close()


def pack(result: dict, files: dict) -> dict:
    """Portable project output snapshot; preserve exact bytes, omit preview copies."""
    stored = copy.deepcopy(result)
    for graphic in stored.get("graphics", []):
        graphic.pop("src", None)
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for key, item in files.items():
            archive.writestr(key, item["data"])
    return {"result": stored, "files": base64.b64encode(buffer.getvalue()).decode("ascii")}


def unpack(snapshot: dict) -> tuple[dict, dict]:
    result = copy.deepcopy(snapshot["result"])
    with zipfile.ZipFile(io.BytesIO(base64.b64decode(snapshot["files"], validate=True))) as archive:
        # Archives are read in memory, never extracted to arbitrary paths.
        if sum(item.file_size for item in archive.infolist()) > 1024 ** 3:
            raise ValueError("Project outputs exceed the 1 GB safety limit")
        files = {info["id"]: {**info, "data": archive.read(info["id"])} for info in result["downloads"]}
    for graphic in result.get("graphics", []):
        item = files[graphic["downloadId"]]
        mime = "image/svg+xml" if item["filename"].endswith(".svg") else "image/png"
        graphic["src"] = f"data:{mime};base64," + base64.b64encode(item["data"]).decode("ascii")
    return result, files


def save(root: Path, result: dict, files: dict) -> dict:
    run = {"id": uuid.uuid4().hex, "created": datetime.now(timezone.utc).isoformat(),
           "title": result.get("project", "Golden Gate design"), "plasmidCount": result.get("plasmidCount", 0)}
    snapshot = pack({**result, "savedRun": run, "temporary": False}, files)
    with _database(root) as db:
        db.execute("INSERT INTO runs (id, created, title, plasmid_count, result_json, files_zip, project_id) VALUES (?, ?, ?, ?, ?, ?, ?)",
                   (run["id"], run["created"], run["title"], run["plasmidCount"], json.dumps(snapshot["result"]), base64.b64decode(snapshot["files"]),
                    (result.get("projectState") or {}).get("projectId")))
    return run


def list_runs(root: Path) -> list[dict]:
    if not (root / "DesignLibrary" / "designs.sqlite3").is_file():
        return []
    with _database(root) as db:
        projects, seen = [], set()
        for row in db.execute("SELECT id, created, title, plasmid_count, project_id FROM runs ORDER BY created DESC, id DESC"):
            key = row[4] or row[0]  # Legacy output-only runs remain accessible.
            if key not in seen:
                projects.append({"id": row[0], "created": row[1], "title": row[2], "plasmidCount": row[3]})
                seen.add(key)
        return projects


def load(root: Path, run_id: str) -> tuple[dict, dict]:
    if not (root / "DesignLibrary" / "designs.sqlite3").is_file():
        raise ValueError("No saved designs yet")
    with _database(root) as db:
        row = db.execute("SELECT result_json, files_zip FROM runs WHERE id = ?", (run_id,)).fetchone()
    if row is None:
        raise ValueError("Saved design not found")
    return unpack({"result": json.loads(row[0]), "files": base64.b64encode(row[1]).decode("ascii")})
