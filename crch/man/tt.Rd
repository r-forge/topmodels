\name{tt}
\alias{dtt}
\alias{ptt}
\alias{rtt}
\alias{qtt}

\title{The Truncated Student-t Distribution}

\description{
  Density, distribution function, quantile function, and random generation
  for the left and/or right truncated student-t distribution with \code{df}
  degrees of freedom.  
}
\usage{
dtt(x, mean = 0, sd = 1, df, left = -Inf, right = Inf, log = FALSE)

ptt(q, mean = 0, sd = 1, df, lower.tail = TRUE, log.p = FALSE, left = -Inf,
  right = Inf)

rtt(n, mean = 0, sd = 1, df, left = -Inf, right = Inf)

qtt(p, mean = 0, sd = 1, df, lower.tail = TRUE, log.p = FALSE, left = -Inf, 
  right = Inf)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
  taken to be the number required.}
  \item{mean}{location parameter.}
  \item{sd}{scale parameter.}
  \item{df}{degrees of freedom (> 0, maybe non-integer).  \code{df = Inf} is
  allowed.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x]
  otherwise, P[X > x].}
  \item{left}{left censoring point.}
  \item{right}{right censoring point.}
}

\details{
If \code{mean} or \code{sd} are not specified they assume the default values
of \code{0} and \code{1}, respectively. \code{left} and \code{right} have the defaults \code{-Inf} and \code{Inf} respectively.

The truncated student-t distribution has density 
\deqn{f(x) = 1/\sigma \tau((x - \mu)/\sigma) / 
  (T((right - \mu)/\sigma) - T((left - \mu)/\sigma))}
for \eqn{left \le x \le right}, and 0 otherwise.
     
where \eqn{T} and \eqn{\tau} are the cumulative distribution function
and probability density function of the student-t distribution with
\code{df} degrees of freedom respectively, \eqn{\mu} is the location of the
distribution, and \eqn{\sigma} the scale.
}

\value{
 \code{dtt} gives the density, \code{ptt} gives the distribution
     function, \code{qtt} gives the quantile function, and \code{rtt}
     generates random deviates.
}

\seealso{\code{\link{dt}}}

\keyword{distribution}
