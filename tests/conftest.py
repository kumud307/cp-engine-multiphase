import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cp_engine import CpEngine  # noqa: E402
from cp_engine.database import load_database  # noqa: E402


@pytest.fixture(scope="session")
def engine():
    return CpEngine()


@pytest.fixture(scope="session")
def db():
    return load_database()
