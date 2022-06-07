


zip:
	-rm abstract_IWSM2022_Stauffer_final.zip
	(cp abstract.tex Stauffer.tex && \
		zip abstract_IWSM2022_Stauffer_final.zip Stauffer.tex Stauffer-rootograms.pdf Stauffer-qqresiduals.pdf)
	rm Stauffer.tex


abstract:
	pdflatex abstract.tex
	pdflatex abstract.tex


.PHONY: clean
clean:
	-rm abstract.aux
	-rm abstract.bbl
	-rm abstract.blg
	-rm abstract.log

.PHONY: distclean
distclean:
	make clean
	-rm abstract.pdf
