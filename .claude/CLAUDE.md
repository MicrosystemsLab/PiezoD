# Overview

Engineering tools for internal use.

## Repository layout

TODO

## Tech Stack

- Python
	- Environment: `uv`
	- Lint: `ruff`
	- Test: `pytest`
	- Type checking: `ty`
	- Build: `hatchling`
- Matlab
	- Lint: `checkcode`
	- Test: `runtests`
- KLayout
	- pcell: use python not ruby

## Commands

- Python
	- Install: `uv sync` from pyproject.toml
	- Lint: `uvx ruff check [path] --fix --extend-select I,B,SIM,C4,ISC,PIE && uvx ruff format [path]`
	- Test: `uv run pytest`
	- Type check: `uvx ty check [path]`
- Matlab
	- Test: `matlab -batch "runtests('tests')"

## Safety

- Never run `git checkout` (destructive)
- Never use unicode in code (breaks Python on Windows)

## Workflow

- Plan first, share plan, get explicit approval before implementing
- Never assume a good plan means you should implement it
- If unsure or stuck, ask the user
- If something fails, do not silently move on - ask for clarification
- If `uv` fails due to permissions, consult user

## Testing

- Run tests after changes; all tests must pass before task is complete
- Prefer tests over manual verification

## Code Quality

- Use type hints and docstrings for all functions
- Keep files under 300 LOC; split if longer
- Fix all linter issues

## Parameters

- Single source of truth for defaults (e.g. config.py)
- Use enums instead of bare strings
- Never use `params.get('key', default)` pattern
- Never use `if config then value else default` pattern

## Refactoring

- Delete old interfaces; no legacy wrappers or thin compatibility layers

## Documentation

- No emojis in code or docs
- No excessive bold in markdown; use styling selectively
- No "last updated" dates or authorship
- Commands on single line, not split across lines

## Hardware

- Read project-specific docs and datasheets before making suggestions
- Use exact part numbers, not generic equivalents
- Prioritize safety over convenience; err on side of caution

## Plotting

- For python, include `addcopyfighandler` for easy matplotlib copy/paste
