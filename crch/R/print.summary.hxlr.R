print.summary.hxlr <-
function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    print(x$info)
    if(NROW(x$coefficients$predictor) | NROW(x$coefficients$intercept)) {
      cat(paste("\nCoefficients:\n", sep = ""))
      printCoefmat(rbind(x$coefficients$intercept, x$coefficients$predictor), digits = digits, signif.legend = FALSE)
    } 

    if(NROW(x$coefficients$scale)) {
      cat(paste("\nlog-scale coefficients:\n", sep = ""))
      printCoefmat(x$coefficients$scale, digits = digits, signif.legend = FALSE)
    } 

    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[, 4L] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
  }
  invisible(x)
}
