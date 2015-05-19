\name{tlogis}
\alias{dtlogis}
\alias{ptlogis}
\alias{rtlogis}
\alias{qtlogis}

\title{The Truncated Logistic Distribution}

\description{
  Density, distribution function, quantile function, and random generation
  for the left and/or right truncated logistic distribution.  
}
\usage{
dtlogis(x, mean = 0, sd = 1, left = -Inf, right = Inf, log = FALSE)

ptlogis(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, left = -Inf,
  right = Inf)

rtlogis(n, mean = 0, sd = 1, left = -Inf, right = Inf)

qtlogis(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE, left = -Inf, 
  right = Inf)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
  taken to be the number required.}
  \item{mean}{location parameter.}
  \item{sd}{scale parameter.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x]
  otherwise, P[X > x].}
  \item{left}{left truncation point.}
  \item{right}{right truncation point.}
}

\details{
If \code{mean} or \code{sd} are not specified they assume the default values
of \code{0} and \code{1}, respectively. \code{left} and \code{right} have the defaults \code{-Inf} and \code{Inf} respectively.

The truncated logistic distribution has density 

\deqn{f(x) = 1/\sigma \lambda((x - \mu)/\sigma) / (\Lambda((right - \mu)/\sigma) - \Lambda((left - \mu)/\sigma))}
for \eqn{left \le x \le right}, and 0 otherwise.
     
\eqn{\Lambda} and \eqn{\lambda} are the cumulative distribution function
and probability density function of the standard logistic distribution
respectively, \eqn{\mu} is the location of the distribution, and \eqn{\sigma} 
the scale.
}

\value{
 \code{dtlogis} gives the density, \code{ptlogis} gives the distribution
     function, \code{qtlogis} gives the quantile function, and \code{rtlogis}
     generates random deviates.
}

\seealso{\code{\link{dlogis}}}

\keyword{distribution}
