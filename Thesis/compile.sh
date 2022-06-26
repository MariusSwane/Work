#!/bin/bash

pdflatex main
bibtex main
bibtex Kappa
bibtex StatePresenceCivilConflict
pdflatex main
pdflatex main
