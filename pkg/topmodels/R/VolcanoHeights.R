#' Tukey's Volcano Heights
#' 
#' Heights of 218 volcanos taken from Tukey (1972).
#' 
#' The data are taken from Tukey (1972) who obtained them from
#' \emph{The World Almanac, 1966} (New York: The New York
#' World-Telegram and The Sun, 1966), pp. 282--283.
#' 
#' @usage data("VolcanoHeights", package = "topmodels")
#' 
#' @format A numeric vector of 218 volcano heights (in 1000 feet).
#' 
#' @source Figure 1 in Tukey (1972).
#' 
#' @references Tukey JW (1972). \dQuote{Some Graphic and Semigraphic Displays.}
#'   In Bancroft TA (ed.), \emph{Statistical Papers in Honor of George W. Snedecor},
#'   pp. 293--316. Iowa State University Press, Ames, IA.
#'   Reprinted in Cleveland WS (ed.): \emph{The Collected Works of John W. Tukey,
#'   Volume V. Graphics: 1965--1985}, Wadsworth & Brooks/Cole, Pacific Grove, CA, 1988.
#'
#' @examples
#' ## Rootograms from Tukey (1972)
#' ## (some 'breaks' don't match exactly)
#' library("topmodels")
#' data("VolcanoHeights", package = "topmodels")
#' 
#' ## Figure 16
#' rootogram(lm(VolcanoHeights ~ 1), style = "standing",
#'   breaks = 0:20 - 0.01, expected = FALSE, confint = FALSE)
#' 
#' ## Figure 17
#' rootogram(lm(sqrt(1000 * VolcanoHeights) ~ 1), style = "standing",
#'   breaks = 0:17 * 10 - 1.1, expected = FALSE, confint = FALSE)
#' 
#' ## Figure 18
#' rootogram(lm(sqrt(1000 * VolcanoHeights) ~ 1), style = "hanging",
#'   breaks = -2:18 * 10 - 1.1, confint = FALSE)
#' 
#' ## Figure 19
#' rootogram(lm(sqrt(1000 * VolcanoHeights) ~ 1), style = "suspended",
#'   breaks = -2:18 * 10 - 1.1, ylim = c(6, -2), confint = FALSE)
#' abline(h = c(-1.5, -1, 1, 1.5), lty = c(2, 3, 3, 2))
"VolcanoHeights"
