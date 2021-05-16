# run `make vignettes && make clean` before `R CMD build ...`

vignettes: vignettes/article.tex vignettes/article2.tex

vignettes/article-original.Rnw: article.Rnw
	cp article.Rnw vignettes/article-original.Rnw

vignettes/article-Sweave.Rtex: vignettes/article-original.Rnw
	R -e 'knitr::knit("vignettes/article-original.Rnw")'
	mv article-original.tex vignettes/article-original.tex
	mv vignettes/article-original.tex vignettes/article-Sweave.Rtex
	mkdir -p vignettes/Figures
	mv -t vignettes/Figures Figures/*.png

# also needed in extrafig: predlik1.pdf predlik3.pdf
# also needed in Figures: jss2014-betapost-1.pdf jss2014-homo-1.pdf jss2014-scatter-1.pdf jss2014-usd2-1.pdf jss2014-usd4-1.pdf jss2014-usd6-1.pdf jss2014-usd8-1.pdf jss2014-hetero-1.pdf jss2014-qqplot-1.pdf jss2014-usd1-1.pdf jss2014-usd3-1.pdf jss2014-usd5-1.pdf jss2014-usd7-1.pdf
vignettes/article.tex: vignettes/article-knitr.Rtex
	cp vignettes/article-knitr.Rtex vignettes/article.tex

vignettes/article-knitr.Rtex: vignettes/article-Sweave.Rtex
	R -e 'knitr::Sweave2knitr("vignettes/article-Sweave.Rtex")'
	mv vignettes/article-Sweave-knitr.Rtex vignettes/article-knitr.Rtex

vignettes/article2-original.Rnw: article2.Rnw
	cp article2.Rnw vignettes/article2-original.Rnw

vignettes/article2-Sweave.Rtex: vignettes/article2-original.Rnw
	R -e 'knitr::knit("vignettes/article2-original.Rnw")'
	mv article2-original.tex vignettes/article2-original.tex
	mv vignettes/article2-original.tex vignettes/article2-Sweave.Rtex
	mkdir -p vignettes/Figures2
	mv -t vignettes/Figures2 Figures2/*.png

# also needed Figures2: jss3813-corimageplot-1.pdf jss3813-evdiag-1.pdf jss3813-loadplot2-1.pdf jss3813-logdetcovn-1.pdf jss3813-preddist-1.pdf jss3813-svlplot-1.pdf jss3813-varplot2-1.pdf jss3813-comtimeplot-1.pdf jss3813-cortimeplot-1.pdf jss3813-fsvprepdata-1.pdf jss3813-loadplot2-2.pdf jss3813-plotsvlpred-1.pdf jss3813-svlbetaplot-1.pdf jss3813-varplot-1.pdf jss3813-voltimeplot-1.pdf
vignettes/article2.tex: vignettes/article2-knitr.Rtex
	cp vignettes/article2-knitr.Rtex vignettes/article2.tex

vignettes/article2-knitr.Rtex: vignettes/article2-Sweave.Rtex
	R -e 'knitr::Sweave2knitr("vignettes/article2-Sweave.Rtex")'
	mv vignettes/article2-Sweave-knitr.Rtex vignettes/article2-knitr.Rtex

clean:
	rm -fv \
		vignettes/article-original.Rnw \
		vignettes/article-knitr.Rtex \
		vignettes/article-Sweave.Rtex
		vignettes/article2-original.Rnw \
		vignettes/article2-knitr.Rtex \
		vignettes/article2-Sweave.Rtex

distclean:
	rm -frv \
		vignettes/article-original.Rnw \
		vignettes/article-knitr.Rtex \
		vignettes/article-Sweave.Rtex \
		vignettes/article.tex \
		vignettes/Figures \
		vignettes/article2-original.Rnw \
		vignettes/article2-knitr.Rtex \
		vignettes/article2-Sweave.Rtex \
		vignettes/article2.tex \
		vignettes/Figures2

