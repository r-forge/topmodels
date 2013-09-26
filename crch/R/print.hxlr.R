print.hxlr <-
function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    if(length(x$coefficients$predictor) | length(x$coefficients$intercept)) {
      cat(paste("Coefficients:\n", sep = ""))
       print.default(format(c(x$coefficients$intercept, x$coefficients$predictor), digits = digits), print.gap = 2, quote = FALSE)
       cat("\n")
    } else cat("No coefficients\n\n")
    if(length(x$coefficients$scale)) {
      cat(paste("Coefficients (scale model with log link):\n", sep = ""))
      print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else cat("No coefficients (in scale model)\n\n")
    cat("---\n")
    cat("Converted coefficients for location and scale:\n")
    if(length(x$coefficients2$location)) {
      cat(paste("Coefficients (location model):\n", sep = ""))
       print.default(format(x$coefficients2$location, digits = digits), print.gap = 2, quote = FALSE)
       cat("\n")
    } else cat("No coefficients (in location model)\n\n")
    if(length(x$coefficients2$scale)) {
      cat(paste("Coefficients (scale model with log link):\n", sep = ""))
      print.default(format(x$coefficients2$scale, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    } else cat("No coefficients (in scale model)\n\n")
    cat("\n")
  }
  invisible(x)
}
