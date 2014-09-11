# type "make" command in Unix to create mets.pdf file 
all:
	pdflatex twoPhaseModeling
	bibtex twoPhaseModeling
	pdflatex twoPhaseModeling
	pdflatex twoPhaseModeling

clean:
	(rm -rf *.spl *.out *.ps *.log *.dvi *.aux *.*% *.lof *.lop *.lot *.toc *.idx *.ilg *.ind *.bbl *blg)
