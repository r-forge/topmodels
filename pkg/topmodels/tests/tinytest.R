# --------------------------------------------------------------------
# TRIGGER TINY TESTS
# --------------------------------------------------------------------
  
if (requireNamespace("tinytest", quietly = TRUE)) {
  tinytest::test_package("topmodels")
}

