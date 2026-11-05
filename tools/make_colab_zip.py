"""Create cp_engine_project.zip for uploading to Google Colab.

Run from the project root:  python tools/make_colab_zip.py
"""
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INCLUDE = ["cp_engine", "data", "validation", "tests", "audit", "notebooks", "docs", "tools", "manuscript", "figures", "pyproject.toml",
           "requirements.txt", "README.md", "MODEL_VERIFICATION_REPORT.md"]
SKIP = {"__pycache__", ".ipynb_checkpoints"}

out = ROOT / "cp_engine_project.zip"
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
    for item in INCLUDE:
        p = ROOT / item
        files = [p] if p.is_file() else [f for f in p.rglob("*") if f.is_file()]
        for f in files:
            if not any(part in SKIP for part in f.parts):
                z.write(f, f.relative_to(ROOT))
print("written", out)
