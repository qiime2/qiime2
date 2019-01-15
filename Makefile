.PHONY: all lint test install dev clean distclean

PYTHON ?= python
PREFIX ?= $(CONDA_PREFIX)

all: ;

lint:
	q2lint
	flake8

test: all
	QIIMETEST= nosetests

install: all
	$(PYTHON) setup.py install && \
	mkdir -p $(PREFIX)/etc/conda/activate.d && \
	cp hooks/00_activate_qiime2_envs.sh $(PREFIX)/etc/conda/activate.d/ && \
	mkdir -p $(PREFIX)/etc/conda/deactivate.d && \
	cp hooks/00_deactivate_qiime2_envs.sh $(PREFIX)/etc/conda/deactivate.d/

dev: all
	pip install -e .

clean: distclean

distclean: ;
