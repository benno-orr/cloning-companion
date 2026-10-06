"""Publish a successfully signed local build to the in-app updater."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mac_app.local_updates import publish, update_directory

if __name__ == "__main__":
    result = publish(Path(sys.argv[1]), update_directory())
    print(f"Published local in-app update {result['version']}")
