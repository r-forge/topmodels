\name{rainIbk}
\alias{rainIbk}
\docType{data}
\title{
Precipitation Observations and Forecasts for Innsbruck
}
\description{
Accumulated 5-8 days precipitation amount for Innsbruck. Data includes 
GEFS reforecasts (Hamill et al. 2012) and observations from SYNOP 
station Innsbruck Airport (11120) from 2000-01-01 to 2013-09-17.
}
\usage{data(rainIbk)}
\format{
  A data frame with 4977 observations on the following 12 variables.
  \describe{
    \item{\code{rain}}{3 days accumulated precipitation amount 
      observations from SYNOP station 11120 (Innsbruck Airport).}
    \item{\code{rainfc.1 ... rainfc.11}}{5-8 days accumulated 
      precipitation amount forecast for SYNOP station 11120 from GEFS 
      reforecast. Each column contains forecasts from one ensemble 
      member.}
  }
}
%\details{
%}
\source{
Observations:
http://www.ogimet.com/synops.phtml.en

Reforecasts: http://www.esrl.noaa.gov/psd/forecasts/reforecast2/
}
\references{
Hamill, T. M., G. T. Bates, J. S. Whitaker, D. R. Murray, M. Fiorino, 
T. J. Galarneau, Jr., Y. Zhu, and W. Lapenta, 2013: NOAA's 
second-generation global medium-range ensemble reforecast data set. 
Bull Amer. Meteor. Soc., in press.
}
\examples{
## Spread skill relationship
## load and prepare data
data(rainIbk)
## mean and standard deviation of square root transformed ensemble 
## forecasts
rainIbk$sqrtensmean <- 
  apply(sqrt(rainIbk[,grep('^rainfc',names(rainIbk))]), 1, mean)
rainIbk$sqrtenssd <- 
  apply(sqrt(rainIbk[,grep('^rainfc',names(rainIbk))]),  1, sd)
  
sdcat <- cut(rainIbk$sqrtenssd, c(-Inf, quantile(rainIbk$sqrtenssd, 
  seq(0.2,0.8,0.2)), Inf))
levels(sdcat) <- c(1:5)
m <- NULL
for(i in levels(sdcat)) {
  m <- c(m, mean((sqrt(rainIbk$rain)[sdcat == i] 
    - rainIbk$sqrtensmean[sdcat == i])^2, na.rm = TRUE))
}
boxplot((sqrt(rain) - sqrtensmean)^2~sdcat, rainIbk, 
  xlab = "Quintile of ensemble standard deviation", 
  ylab = "mean squared error", main = "Spread skill relationship")
}
\keyword{datasets}
