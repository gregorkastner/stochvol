vignettes: vignettes/article2.Rtex

vignettes/article2-original.Rnw: article2.Rnw
	cp article2.Rnw vignettes/article2-original.Rnw

vignettes/article2-Sweave.Rtex: vignettes/article2-original.Rnw
	R -e 'knitr::knit("vignettes/article2-original.Rnw")'
	mv article2-original.tex vignettes/article2-original.tex
	mv vignettes/article2-original.tex vignettes/article2-Sweave.Rtex

vignettes/article2.Rtex: vignettes/article2-knitr.Rtex
	cp vignettes/article2-knitr.Rtex vignettes/article2.Rtex

vignettes/article2-knitr.Rtex: vignettes/article2-Sweave.Rtex
	R -e 'knitr::Sweave2knitr("vignettes/article2-Sweave.Rtex")'
	mv vignettes/article2-Sweave-knitr.Rtex vignettes/article2-knitr.Rtex

clean:
	rm \
		vignettes/article2-original.Rnw \
		vignettes/article2-knitr.Rtex \
		vignettes/article2-Sweave.Rtex

distclean:
	rm \
		vignettes/article2-original.Rnw \
		vignettes/article2-knitr.Rtex \
		vignettes/article2-Sweave.Rtex \
		vignettes/article2.Rtex

