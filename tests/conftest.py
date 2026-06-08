import sys
from pathlib import Path

import pytest
from click.testing import CliRunner

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT_STR = str(PROJECT_ROOT)
TESTS_ROOT = Path(__file__).resolve().parent

if PROJECT_ROOT_STR not in sys.path:
    sys.path.insert(0, PROJECT_ROOT_STR)


@pytest.fixture(scope="session")
def tests_root() -> Path:
    return TESTS_ROOT


@pytest.fixture(scope="session")
def fixtures_root(tests_root: Path) -> Path:
    return tests_root / "fixtures"


@pytest.fixture(scope="session")
def count_fixtures_root(fixtures_root: Path) -> Path:
    return fixtures_root / "count"


@pytest.fixture
def minimal_count_input(count_fixtures_root: Path) -> Path:
    return count_fixtures_root / "minimal_test_input.tsv"


@pytest.fixture
def cli_runner() -> CliRunner:
    return CliRunner()
