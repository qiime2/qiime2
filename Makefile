.PHONY: all lint test install dev clean distclean

all: ;

lint:
	q2lint
	flake8

test: all
	QIIMETEST= nosetests

install: all
	python setup.py install

dev: all
	pip install -e .

clean: distclean

distclean: ;
