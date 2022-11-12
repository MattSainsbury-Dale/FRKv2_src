#!/bin/bash

mkdir -p ./arXiv_submission

cp jss.cls arXiv_submission/jss.cls
cp jss.bst arXiv_submission/jss.bst
cp FRKv2.tex arXiv_submission/FRKv2.tex
cp FRKv2.bib arXiv_submission/FRKv2.bib
cp FRKv2.bbl arXiv_submission/FRKv2.bbl
cp -r results arXiv_submission/results

## Replace the first line (\documentclass[article]{jss}) with the "nojss" option
## NB: single quotes are important as they allow \ to be escaped
## NB: only works on mac
var='\\documentclass[nojss]{jss}'
sed -i '' "1s/.*/$var/" arXiv_submission/FRKv2.tex

zip -r arXiv_submission ./arXiv_submission
