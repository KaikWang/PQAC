.PHONY: describe check-lazer same-x-demo memver-quadratic-demo acc-memver-layout export-acc-memver import-acc-memver-check import-acc-memver-labrador

PYTHON ?= /home/wkk/miniconda3/envs/sage/bin/python
SDE ?= /home/wkk/tools/sde-external-10.8.0-2026-03-15-lin/sde64 -icx --

describe:
	$(PYTHON) revac_show_relation.py

check-lazer:
	$(PYTHON) revac_show_labrador_skeleton.py --check

same-x-demo:
	$(PYTHON) revac_same_x_demo.py

memver-quadratic-demo:
	$(SDE) $(PYTHON) revac_memver_quadratic_demo.py --primesize 38

acc-memver-layout:
	$(PYTHON) revac_acc_memver_full_layout.py

export-acc-memver:
	../build/revac_acc_memver_export > /tmp/revac_memver.json

import-acc-memver-check: export-acc-memver
	$(PYTHON) revac_acc_memver_import_labrador.py /tmp/revac_memver.json

import-acc-memver-labrador: export-acc-memver
	$(SDE) $(PYTHON) revac_acc_memver_import_labrador.py /tmp/revac_memver.json --mode labrador-simple
