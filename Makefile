SHELL := /bin/bash

PYTHON_BIN ?= python3
LATEXMK ?= latexmk

.DEFAULT_GOAL := check
.PHONY: check check-asymptotic check-finite build build-asymptotic build-finite

check: check-asymptotic check-finite

check-asymptotic:
	"$(PYTHON_BIN)" asymptotic_single_pass_floor/verify_constants.py

check-finite:
	env PYTHON_BIN="$(PYTHON_BIN)" finite_n4096_certificate/code/run_all.sh

build: build-asymptotic build-finite

build-asymptotic:
	cd asymptotic_single_pass_floor && "$(LATEXMK)" -pdf -interaction=nonstopmode -halt-on-error main.tex

build-finite:
	"$(PYTHON_BIN)" finite_n4096_certificate/figures/generate_figures.py
	cd finite_n4096_certificate && "$(LATEXMK)" -pdf -interaction=nonstopmode -halt-on-error main.tex
	mkdir -p finite_n4096_certificate/output/pdf
	cp finite_n4096_certificate/main.pdf finite_n4096_certificate/output/pdf/forr4_n4096_advantage.pdf
