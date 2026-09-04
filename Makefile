VENV ?= .venv
PYTHON ?= $(VENV)/bin/python
BOOTSTRAP_PYTHON ?= python3

.PHONY: venv bootstrap setup dev-install lab test lint typecheck check

venv:
	$(BOOTSTRAP_PYTHON) -m venv $(VENV)

bootstrap: venv
	$(PYTHON) -m pip install --upgrade pip setuptools wheel

setup: bootstrap
	$(PYTHON) -m pip install -e .

dev-install: bootstrap
	$(PYTHON) -m pip install -e .[dev]

lab:
	$(PYTHON) -m jupyter lab

test:
	$(PYTHON) -m pytest

lint:
	$(PYTHON) -m ruff check src tests scripts

typecheck:
	$(PYTHON) -m mypy src tests

check: test lint typecheck
