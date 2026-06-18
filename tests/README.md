# Tests

This directory contains the pytest suite for workflow scripts and supporting fixtures.

## Layout

- `tests/conftest.py` keeps shared pytest fixtures and makes the repo root importable.
- `tests/<area>/test_*.py` holds the actual tests, grouped by script area such as `count`.
- `tests/fixtures/<area>/` stores reusable input data for those tests.

## Adding a new script test

1. Put the new test in the matching area folder, for example `tests/count/test_new_script.py`.
2. Add any reusable inputs under `tests/fixtures/<area>/`.
3. Prefer `click.testing.CliRunner` for Click commands and pytest fixtures like `tmp_path` for temporary outputs.
4. Run a focused check with `conda run -n mpralib python -m pytest -q tests/<area>/test_new_script.py`.

## Current example

The MPRAnalyze compiler test lives in `tests/count/test_mpranalyze_compiler.py` and uses the fixture input at `tests/fixtures/count/minimal_test_input.tsv`.
