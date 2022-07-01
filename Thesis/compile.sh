#!/bin/bash

pdflatex main
bibtex main
bibtex Kappa/Kappa
bibtex ARC/main_text_jpr_final
bibtex BeyondEthnicity/hse_main
bibtex StatePresenceCivilConflict/StatePresenceCivilConflict
bibtex CommunalViolence/OMT
pdflatex main
pdflatex main
