#!/bin/bash

rm MRIbook.aux MRIbook.log MRIbook.toc MRIbook.out 
pdflatex MRIbook.tex
pdflatex MRIbook.tex
pdflatex MRIbook.tex
rm MRIbook.aux MRIbook.log MRIbook.out 
