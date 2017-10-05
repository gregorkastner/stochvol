all: pdf

pdf: article.pdf heavytails.pdf

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

clean:
	rm -rf article.tex article.aux article.bbl article.blg article.log article.out article-*.pdf heavytails.tex heavytails.aux heavytails.bbl heavytails.blg heavytails.log heavytails.out heavytails-*.pdf
