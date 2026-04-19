# Repository Guidelines

## Project Structure & Module Organization
This repository currently contains a single Python virtual environment at `synthetix/` and no checked-in application source, tests, or assets yet. Keep environment files isolated inside `synthetix/`, and add project code at the repository root when development begins.

Recommended layout for future work:
- `src/` for package code
- `tests/` for automated tests
- `docs/` for notes and design references
- `assets/` for static inputs such as sample data or images

Avoid committing files from `synthetix/lib/` or other generated virtualenv contents unless the repository is explicitly meant to vendor dependencies.

## Build, Test, and Development Commands
Use the local virtual environment before running Python tooling:

```bash
source synthetix/bin/activate
python --version
pip list
```

- `source synthetix/bin/activate`: activates the repository’s Python 3.11 environment
- `python --version`: confirms the interpreter in use
- `pip list`: shows installed packages for debugging dependency state

If application code is added later, prefer standard entry points such as `pytest`, `python -m package_name`, and formatter/linter commands documented in `pyproject.toml`.

## Coding Style & Naming Conventions
Use 4-space indentation and follow PEP 8 for Python code. Prefer:
- `snake_case` for modules, files, and functions
- `PascalCase` for classes
- `UPPER_SNAKE_CASE` for constants

Keep modules small and place reusable logic in `src/` instead of notebooks or ad hoc scripts. If formatting tools are introduced, standardize on one source of truth in `pyproject.toml`.

## Testing Guidelines
No test framework is configured yet. When tests are added, use `pytest` with files named `test_*.py` under `tests/`. Keep one test module per source module where practical, and cover both expected behavior and failure cases.

Run tests with:

```bash
source synthetix/bin/activate
pytest
```

## Commit & Pull Request Guidelines
Local Git metadata is not present in this workspace, so no repository-specific commit convention can be inferred. Use short, imperative commit messages such as `Add parser for sample dataset`.

For pull requests, include:
- a clear summary of the change
- setup or migration notes if dependencies changed
- test evidence, such as `pytest` output
- screenshots only when UI or visual assets are introduced
