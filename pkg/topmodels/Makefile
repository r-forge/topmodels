



# ---------------------------------------------------
# Devtools based cmds
# ---------------------------------------------------
document:
	Rscript -e "library('scoringRules'); devtools::document()"

install: document
	Rscript -e "devtools::install()"

check: document
	Rscript -e "library('distributions3'); library('scoringRules'); devtools::check()"

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
