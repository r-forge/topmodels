
library("crch")

uibkblue <- rgb(18, 50, 111, maxColorValue = 255)
uibkorange <- rgb(242, 146, 0, maxColorValue = 255)


load("data.rda")

data <- na.omit(data)


linear <- crch(Tmax~tmax2m, data, dist = "logistic")


## plot functions that are used frequently
scatter <- function() {
  plot(Tmax~tmax2m, data, pch = 19, col = gray(0.2, alpha = 0.1), xlab = "ECMWF forecast [degree C]", ylab = "observation [degree C]", main = "Innsbruck JJA maximum temperature (+30h to +42h)", xlim = c(min(data, na.rm = TRUE), max(data, na.rm = TRUE)), ylim = c(min(data, na.rm = TRUE), max(data, na.rm = TRUE)))
#  abline(h = 0, lwd = 0.5)
#  abline(v = 0, lwd = 0.5)
}

addpdf <- function(point, scalefactor = 1, legend = TRUE, ...) {
  loc <- predict(linear, newdata = data.frame(tmax2m = point))
  scale <- predict(linear, newdata = data.frame(tmax2m = point), type = "scale") * scalefactor
  lines(rep(point, 2), c(loc - 11, loc + 11), lty = 2)
  lines(c(point, point + 5), c(loc, loc), lty = 2)
  lines(dlogis(seq(loc-10, loc+10, 0.1), loc, scale)*20 + point, seq(loc-10, loc+10, 0.1), col = uibkblue, lwd = 2, ...)
  if(legend) legend("bottomright", lwd = 2, col = uibkorange, legend = expression(mu *" = "* beta[0] + beta[1] * bar(ens)), bty = "n")
}

addlik <- function(point, lower = NULL, upper = NULL) {
  loc <- predict(linear, newdata = data.frame(tmax2m = point))
  scale <- predict(linear, newdata = data.frame(tmax2m = point), type = "scale")
  abline(h = upper)
  abline(h = lower)
  if(is.null(lower)) lower <- loc - 10
  if(is.null(upper)) upper <- loc + 10
  x <- c(point, point, dlogis(seq(upper, lower, -0.1), loc, scale)*20 + point)
  y <- c(lower, upper, seq(upper, lower, -0.1))
  polygon(x, y, col = adjustcolor(uibkblue, alpha = 0.7), border = FALSE)
}

addlik2 <- function(point, obs) {
  loc <- predict(linear, newdata = data.frame(tmax2m = point))
  scale <- predict(linear, newdata = data.frame(tmax2m = point), type = "scale")
  lines(c(point, point + dlogis(obs, loc, scale)*20), c(obs, obs), col = uibkblue, lwd = 2)
  points(point, obs, pch = 19)
}

