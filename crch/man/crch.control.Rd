\name{crch.control}
\alias{crch.control}

\title{Auxiliary Function for Controlling crch Fitting}

\description{
  Auxiliary function for \code{crch} fitting. Specifies a list of values passed
  to \code{\link{optim}}.
}
\usage{
crch.control(method = "BFGS", maxit = 5000, hessian = NULL,
  trace = FALSE, start = NULL, dot = "separate", \dots)}
\arguments{
  \item{method}{optimization method passed to \code{\link{optim}}}
  \item{maxit}{the maximum number of iterations.}
  \item{hessian}{logical. Should the numerical Hessian matrix from the 
    \code{\link{optim}} output be used for estimation of the covariance matrix? 
    For \code{dist="gaussian"}, \code{dist="logistic"}, and \code{dist="student"}
    with predefined \code{df} analytical solution is employed by default. 
    For \code{dist="student"} without prespecified \code{df}, no analytical 
    solution is available and a numerical Hessian matrix is forced.}
  \item{trace}{non-negative integer. If positive, tracing information on the
  progress of the optimization is produced.}
  \item{start}{initial values for the parameters to be optimized over.}
  \item{dot}{character specifying how to process formula parts with a dot
    (\code{.}) on the right-hand side. This can either be \code{"separate"}
    so that each formula part is expanded separately or \code{"sequential"}
    so that the parts are expanded sequentially conditional on all prior parts.
    Default is \code{"separate"}}
  \item{\dots}{Additional parameters passed to \code{\link{optim}}. }
}

\value{
  A list with components named as the arguments.
}

\seealso{\code{\link{crch}}, \code{\link{optim}}}

\keyword{regression}
