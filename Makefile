PYTHON ?= python

.PHONY: help install doctor smoke verify figures full test integrity release
help:
	@echo "Targets: install doctor smoke verify figures full test integrity release"
install:
	$(PYTHON) -m pip install -r requirements-lock.txt
doctor:
	$(PYTHON) reproduce.py doctor
smoke:
	$(PYTHON) reproduce.py smoke
verify:
	$(PYTHON) reproduce.py verify
figures:
	$(PYTHON) reproduce.py figures
full:
	$(PYTHON) reproduce.py full
test:
	$(PYTHON) -m unittest discover -s tests -v
integrity:
	$(PYTHON) reproduce.py integrity --release
release:
	$(PYTHON) tools/build_release.py
