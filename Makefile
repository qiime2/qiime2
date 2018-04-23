.PHONY: all lint test install dev clean distclean

all: ;

lint:
	q2lint
	flake8

test: all
	QIIMETEST= nosetests

install: all
	python setup.py install && \
	mkdir -p $(CONDA_PREFIX)/etc/conda/activate.d && \
	cp scripts/activate_qiime2_envs.sh $(CONDA_PREFIX)/etc/conda/activate.d/ && \
	mkdir -p $(CONDA_PREFIX)/etc/conda/deactivate.d && \
	cp scripts/deactivate_qiime2_envs.sh $(CONDA_PREFIX)/etc/conda/deactivate.d/

dev: all
	pip install -e .

clean: distclean

distclean: ;
