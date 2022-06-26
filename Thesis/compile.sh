#!/bin/bash

pdflatex main
bibtex main
bibtex Kappa/Kappa
bibtex StatePresenceCivilConflict/StatePresenceCivilConflict
bibtex CommunalViolence/OMT
pdflatex main
pdflatex main
