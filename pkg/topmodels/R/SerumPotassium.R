#' Serum Potassium Levels
#' 
#' Sample of 152 serum potassium levels.
#' 
#' The data are taken from Rice (2007) who obtained
#' the data from Martin, Gudzinowicz and Fanger (1975) 
#' and reports them rounded to one digit.
#' 
#' @usage data("SerumPotassium", package = "topmodels")
#' 
#' @format A numeric vector of 152 serum potassium levels.
#' 
#' @source Page 350 in Rice (2007).
#' 
#' @references Rice JA (2007). \emph{Mathematical Statistics and Data Analysis},
#'   3rd ed. Duxbury, Belmont, CA.
#' 
#' Martin HF, Gudzinowicz BJ, Fanger H (1975).
#'   \emph{Normal Values in Clinical Chemistry: A Guide to Statistical Analysis of Laboratory Data}.
#'   Marcel Dekker, New York.
#' 
#' @examples
#' library("topmodels")
#' data("SerumPotassium", package = "topmodels")
#' 
#' ## Figure 9.3a-c from Rice (2007), and actual hanging rootogram
#' ## (note that Rice erroneously refers to suspended rootograms as hanging)
#' sp <- lm(SerumPotassium ~ 1)
#' br <- 32:54/10 - 0.05
#' rootogram(sp, scale = "raw", style = "standing",
#'   breaks = br, col = "transparent")
#' rootogram(sp, scale = "raw", style = "suspended",
#'   breaks = br, col = "transparent", ylim = c(2.8, -4))
#' rootogram(sp, scale = "sqrt", style = "suspended",
#'   breaks = br, col = "transparent", ylim = c(1, -1.5))
#' rootogram(sp, breaks = br)
"SerumPotassium"
