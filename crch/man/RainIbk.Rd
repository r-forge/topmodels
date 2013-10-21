\name{RainIbk}
\alias{RainIbk}
\docType{data}
\title{
Precipitation Observations and Forecasts for Innsbruck
}
\description{
Accumulated 5-8 days precipitation amount for Innsbruck. Data includes GEFS
reforecasts (Hamill et al. 2012) and observations from SYNOP station Innsbruck
Airport (11120) from 2000-01-01 to 2013-09-17.
}
\usage{data(RainIbk)}
\format{
  A data frame with 4977 observations on the following 12 variables.
  \describe{
    \item{\code{rain}}{3 days accumulated precipitation amount observations
    from SYNOP station 11120 (Innsbruck Airport).}
    \item{\code{rainfc.1}}{5-8 days accumulated precipitation amount forecast
    for SYNOP station 11120 from GEFS reforecast ensemble member 1.}
    \item{\code{rainfc.2}, \code{rainfc.3}, \code{rainfc.4}, \code{rainfc.5},
    \code{rainfc.6}, \code{rainfc.7}, \code{rainfc.8}, \code{rainfc.9},
    \code{rainfc.10}}{as \code{rainfc.1} but for ensemble members 2-10}
  }
}
\source{
Observations: http://www.ogimet.com/synops.phtml.en

Reforecasts: http://www.esrl.noaa.gov/psd/forecasts/reforecast2/
}
\references{
Hamill, T. M., G. T. Bates, J. S. Whitaker, D. R. Murray, M. Fiorino, 
T. J. Galarneau, Jr., Y. Zhu, and W. Lapenta, 2013: NOAA's 
second-generation global medium-range ensemble reforecast data set. 
Bull Amer. Meteor. Soc., in press.
}
\examples{
## Spread skill relationship ##

## load and prepare data
data(RainIbk)

## mean and standard deviation of square root transformed ensemble 
## forecasts
RainIbk$sqrtensmean <- 
  apply(sqrt(RainIbk[,grep('^rainfc',names(RainIbk))]), 1, mean)
RainIbk$sqrtenssd <- 
  apply(sqrt(RainIbk[,grep('^rainfc',names(RainIbk))]),  1, sd)

## quintiles of sqrtenssd
sdcat <- cut(RainIbk$sqrtenssd, c(-Inf, quantile(RainIbk$sqrtenssd, 
  seq(0.2,0.8,0.2)), Inf), labels = c(1:5))

## mean forecast errors for each quintile
m <- NULL
for(i in levels(sdcat)) {
  m <- c(m, mean((sqrt(RainIbk$rain)[sdcat == i] 
    - RainIbk$sqrtensmean[sdcat == i])^2, na.rm = TRUE))
}

## plot
boxplot((sqrt(rain) - sqrtensmean)^2~sdcat, RainIbk, 
  xlab = "Quintile of ensemble standard deviation", 
  ylab = "mean squared error", main = "Spread skill relationship")

}
\keyword{datasets}
