.PHONY: all lint test install dev clean distclean

PYTHON ?= python
PREFIX ?= $(CONDA_PREFIX)

all: ;

lint:
	q2lint
	flake8

test: all
	QIIMETEST= pytest --doctest-modules

jupyter-ext: all
	mkdir -p $(PREFIX)/etc/jupyter/jupyter_server_config.d/ && \
	cp qiime2/jupyter/qiime2.json $(PREFIX)/etc/jupyter/jupyter_server_config.d/

# for parallel, pip install pytest-xdist
mystery-stew: all
	MYSTERY_STEW= pytest qiime2/tests/mystery_stew.py -n auto

install: all jupyter-ext
	$(PYTHON) -m pip install -v . && \
	mkdir -p $(PREFIX)/etc/conda/activate.d && \
	cp hooks/00_activate_qiime2_envs.sh $(PREFIX)/etc/conda/activate.d/ && \
	mkdir -p $(PREFIX)/etc/conda/deactivate.d && \
	cp hooks/00_deactivate_qiime2_envs.sh $(PREFIX)/etc/conda/deactivate.d/

dev: all jupyter-ext
	pip install -e .

clean: distclean

distclean: ;
