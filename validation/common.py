"""Shared helpers for the validation scripts.

AARD (average absolute relative deviation) and maximum deviation:

    AARD % = (100/N) Σ |y_model − y_exp| / |y_exp|
    max  % = 100 max |y_model − y_exp| / |y_exp|
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

DATA = Path(__file__).resolve().parent / "data"
RESULTS = Path(__file__).resolve().parent / "results"
FIGS = ROOT / "figures"


def aard(model, exp) -> tuple[float, float, int]:
    """Return (AARD %, max |dev| %, N)."""
    model, exp = np.asarray(model, dtype=float), np.asarray(exp, dtype=float)
    dev = np.abs(model - exp) / np.abs(exp) * 100.0
    return float(dev.mean()), float(dev.max()), int(dev.size)


def load(name: str) -> pd.DataFrame:
    """Read a CSV in validation/data (lines starting with '#' are the SOURCE/TODO header)."""
    df = pd.read_csv(DATA / name, comment="#")
    return df.dropna(how="all")


def row(case: str, system: str, quantity: str, n: int, mean: float | None, maximum: float | None,
        reference: str, status: str = "ok", metric: str = "AARD %") -> dict:
    """One line of the validation summary. `metric` says what mean/max are:
    'AARD %' (relative deviation), 'RMSD J/mol', '|ΔT| K' or 'rel. dev %'."""
    return {"case": case, "system": system, "quantity": quantity, "metric": metric, "N": n,
            "mean": mean, "max": maximum, "reference": reference, "status": status}


def pending(case: str, system: str, quantity: str, csv: str) -> dict:
    msg = f"no data yet - paste literature data into validation/data/{csv}"
    print(f"  [{case}] {system}: {msg}")
    return row(case, system, quantity, 0, None, None, "see docs/DATA_TO_VERIFY.md", "PENDING DATA")
