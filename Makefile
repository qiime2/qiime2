.PHONY: all lint test clean

all: ;

lint:
	q2lint
	flake8

test: all
	nosetests

test-cov: all
	nosetests

install: all
	python setup.py install

dev: all
	pip install -e .

clean: distclean

distclean: ;
