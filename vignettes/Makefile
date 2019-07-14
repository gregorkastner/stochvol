all: pdf

pdf: article.pdf heavytails.pdf paper.pdf

article.pdf: article.tex
	pdflatex article
	bibtex article
	pdflatex article
	pdflatex article

heavytails.pdf: heavytails.tex
	pdflatex heavytails
	bibtex heavytails
	pdflatex heavytails
	pdflatex heavytails

paper.pdf: paper.tex
	pdflatex paper
	bibtex paper
	pdflatex paper
	pdflatex paper

clean:
	rm -rfv ./*.{aux,bbl,blg,fdb_latexmk,fls,log,out,tex}

distclean: clean
	rm -rfv Figures
	rm -fv paper.Rnw
	rm -fv ref.bib
