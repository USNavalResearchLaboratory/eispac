
MAIN = eag

all:
	pdflatex $(MAIN) 
	bibtex $(MAIN) 
	pdflatex $(MAIN) 
	pdflatex $(MAIN) 
	cp $(MAIN).pdf EISPAC_Users_Guide.pdf
	open EISPAC_Users_Guide.pdf

clean:
	rm -rf *.aux *.log *.bbl *.blg *.dvi *.pdf *.toc

