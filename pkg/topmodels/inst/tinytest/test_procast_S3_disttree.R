# --------------------------------------------------------------------
# TESTS FOR FUNCTIONS WITHIN `procast.R` FOR S3-CLASS `crch`
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Test procast.crch argument = `type = "parameter"`
# --------------------------------------------------------------------
require("disttree")
dt1 <- disttree(dist ~ speed, data = cars)
dt2 <- disttree(dist ~ speed, data = cars, family = dist_list_normal)

df1 <- distforest(dist ~ speed, data = cars)
df2 <- distforest(dist ~ speed, data = cars, family = dist_list_normal)

expect_equal(
  procast(dt1, type = "parameter"),
  predict(dt1, type = "parameter"),
  check.attributes = FALSE
)

expect_equal(
  procast(dt2, type = "parameter"),
  predict(dt2, type = "parameter"),
  check.attributes = FALSE
)

expect_equal(
  procast(df1, type = "parameter"),
  predict(df1, type = "parameter"),
  check.attributes = FALSE
)

expect_equal(
  procast(df2, type = "parameter"),
  predict(df2, type = "parameter"),
  check.attributes = FALSE
)


# --------------------------------------------------------------------
# Test for NO() consistency of <p,d,q>norm() with distfamily: Tree
# --------------------------------------------------------------------
expect_equal(
  procast(dt1, at = c(0.25, 0.5, 0.75), type = "quantile"),
  procast(dt1, at = c(0.25, 0.5, 0.75), type = "quantile", use_distfamily = FALSE),
)

expect_equal(
  procast(dt1, at = c(0.25, 0.5, 0.75), type = "location"),
  procast(dt1, at = c(0.25, 0.5, 0.75), type = "location", use_distfamily = FALSE),
)

expect_equal(
  procast(dt1, at = c(0.25, 0.5, 0.75), type = "scale"),
  procast(dt1, at = c(0.25, 0.5, 0.75), type = "scale", use_distfamily = FALSE),
)

expect_equal(
  procast(dt1, at = c(0.25, 0.5, 0.75), type = "parameter"),
  procast(dt1, at = c(0.25, 0.5, 0.75), type = "parameter", use_distfamily = FALSE),
)

expect_equal(
  procast(dt1, at = c(0.25, 0.5, 0.75), type = "density", log = FALSE),
  procast(dt1, at = c(0.25, 0.5, 0.75), type = "density", use_distfamily = FALSE),
)

expect_equal(
  procast(dt1, at = c(0.25, 0.5, 0.75), type = "probability"),
  procast(dt1, at = c(0.25, 0.5, 0.75), type = "probability", use_distfamily = FALSE),
)


# --------------------------------------------------------------------
# Test for NO() consistency of <p,d,q>norm() with distfamily: Forest
# --------------------------------------------------------------------
expect_equal(
  procast(df1, at = c(0.25, 0.5, 0.75), type = "quantile"),
  procast(df1, at = c(0.25, 0.5, 0.75), type = "quantile", use_distfamily = FALSE),
)

expect_equal(
  procast(df1, at = c(0.25, 0.5, 0.75), type = "location"),
  procast(df1, at = c(0.25, 0.5, 0.75), type = "location", use_distfamily = FALSE),
)

expect_equal(
  procast(df1, at = c(0.25, 0.5, 0.75), type = "scale"),
  procast(df1, at = c(0.25, 0.5, 0.75), type = "scale", use_distfamily = FALSE),
)

expect_equal(
  procast(df1, at = c(0.25, 0.5, 0.75), type = "parameter"),
  procast(df1, at = c(0.25, 0.5, 0.75), type = "parameter", use_distfamily = FALSE),
)

expect_equal(
  procast(df1, at = c(0.25, 0.5, 0.75), type = "density", log = FALSE),
  procast(df1, at = c(0.25, 0.5, 0.75), type = "density", use_distfamily = FALSE),
)

expect_equal(
  procast(df1, at = c(0.25, 0.5, 0.75), type = "probability"),
  procast(df1, at = c(0.25, 0.5, 0.75), type = "probability", use_distfamily = FALSE),
)


# --------------------------------------------------------------------
# Test for NO() consistency of <p,d,q>norm() with distfamily: Tree w/ newdata
# --------------------------------------------------------------------
nd <- data.frame(speed = c(10, 15, 20))

expect_equal(
  procast(dt1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "quantile"),
  procast(dt1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "quantile", use_distfamily = FALSE),
)

expect_equal(
  procast(dt1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "location"),
  procast(dt1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "location", use_distfamily = FALSE),
)

expect_equal(
  procast(dt1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "scale"),
  procast(dt1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "scale", use_distfamily = FALSE),
)

expect_equal(
  procast(dt1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "parameter"),
  procast(dt1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "parameter", use_distfamily = FALSE),
)

expect_equal(
  procast(dt1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "density", log = FALSE),
  procast(dt1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "density", use_distfamily = FALSE),
)

expect_equal(
  procast(dt1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "probability"),
  procast(dt1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "probability", use_distfamily = FALSE),
)

# --------------------------------------------------------------------
# Test for NO() consistency of <p,d,q>norm() with distfamily: Forest w/ newdata
# --------------------------------------------------------------------
nd <- data.frame(speed = c(10, 15, 20))

expect_equal(
  procast(df1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "quantile"),
  procast(df1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "quantile", use_distfamily = FALSE),
)

expect_equal(
  procast(df1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "location"),
  procast(df1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "location", use_distfamily = FALSE),
)

expect_equal(
  procast(df1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "scale"),
  procast(df1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "scale", use_distfamily = FALSE),
)

expect_equal(
  procast(df1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "parameter"),
  procast(df1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "parameter", use_distfamily = FALSE),
)

expect_equal(
  procast(df1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "density", log = FALSE),
  procast(df1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "density", use_distfamily = FALSE),
)

expect_equal(
  procast(df1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "probability"),
  procast(df1, at = c(0.25, 0.5, 0.75), newdata = nd, type = "probability", use_distfamily = FALSE),
)
