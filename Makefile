VENV ?= .venv-dev
PYTHON ?= $(VENV)/bin/python
BOOTSTRAP_PYTHON ?= python3.12
CONSTRAINTS ?= constraints-dev.txt

.PHONY: venv bootstrap setup dev-install lab test lint typecheck check

venv:
	$(BOOTSTRAP_PYTHON) -m venv $(VENV)

bootstrap: venv
	$(PYTHON) -m pip install --upgrade pip setuptools wheel

setup: bootstrap
	$(PYTHON) -m pip install -e .

dev-install: venv
	$(PYTHON) -c 'import sys; assert sys.version_info[:2] == (3, 12), "Development and CI require Python 3.12"'
	$(PYTHON) -m pip install -c $(CONSTRAINTS) pip setuptools wheel
	$(PYTHON) -m pip install --no-build-isolation -c $(CONSTRAINTS) -e '.[dev]'
	$(PYTHON) -m pip check

lab:
	$(PYTHON) -m jupyter lab

test:
	$(PYTHON) -m pytest

lint:
	$(PYTHON) -m ruff check src tests scripts

typecheck:
	$(PYTHON) -m mypy src tests

check: test lint typecheck
