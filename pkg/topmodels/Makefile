


# Extracting current package version
VERSION := $(shell grep '^Version:' DESCRIPTION | awk '{print $$2}')

# ---------------------------------------------------
# Devtools based cmds
# ---------------------------------------------------
document:
	Rscript -e "library('scoringRules'); devtools::document()"

.PHONY: install build clean
.PHONY: clean
clean:
	-rm src/*.so
	-rm src/*.o
	-rm ../topmodels_*.tar.gz
build: clean document
	@echo Building current version: $(VERSION)
	(cd ../ && R CMD build topmodels)
install: build
	@echo Installing current version: $(VERSION)
	(cd ../ && R CMD INSTALL topmodels_$(VERSION).tar.gz)
check: build
	@echo Checking current version: $(VERSION)
	(cd ../ && R CMD check --as-cran topmodels_$(VERSION).tar.gz)


test: install
	Rscript -e "tinytest::test_package('topmodels')"
testwarn: install
	Rscript -e "tinytest::test_package('topmodels'); warnings()"

bamlsstest: install
	Rscript -e "tinytest::test_package('topmodels', testdir = 'tinytest_bamlss')"

coverage: install
	Rscript -e 'covr::report(covr::package_coverage(), file = "../coverage.html")'

reto:
	#Rscript -e "library('topmodels'); devtools::load_all('topmodels'); tinytest::run_test_file('topmodels/inst/tinytest/test_pithist_usage.R')"
	Rscript -e "library('topmodels'); devtools::load_all(); tinytest::run_test_file('inst/tinytest/test_pithist_02_values.R')"


# ---------------------------------------------------
# Command line version
# ---------------------------------------------------
packageversion:=$(shell cat DESCRIPTION | egrep Version | sed 's/Version://g')
cmd: SHELL:=/bin/bash
cmd: cmd
	(cd ../ && \
		R CMD build --no-build-vignettes topmodels && \
		R CMD INSTALL topmodels_$(shell printf "%s"${packageversion}).tar.gz)
