\name{predict.hxlr}
\alias{predict.hxlr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  predict hxlr models
}
\description{
Predict function for heteroscedastic extended logistic regression models
}
\usage{
\method{predict}{hxlr}(object, newdata = object$data, type = c("class", "probability", "cumprob", "location", "scale"), thresholds = object$thresholds, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    fitted object of class 'hxlr'
}
  \item{newdata}{
   An optional data frame in which to look for variables with which to predict.  If omitted, same data as for fitting is used.
}
  \item{type}{
    Type of prediction: '"prob"' returns data frame with category probabilities, '"cumprob"' returns cumulative probabilities, '"location"' and "'scale"' return the location and scale of the predicted latent distribution respectively, and "'class"' returns the category with the highest probability. Default is "'class"'.
}
  \item{thresholds}{
    Optional thresholds used for defining the thresholds (for types '"prob"', '"cumprob"', and '"class"'). Can be different from thresholds used for fitting. If omitted, the same thresholds as for fitting are used.
}
  \item{...}{
    Further arguments passed to or from other methods.
}
}
%\details{
%TODO
%}
\value{
For type '"prob"' a matrix with number of intervals (= number of thresholds + 1) columns is produced. Each row corresponds to a row in newdata and contains the predicted probabilities to fall in the corresponding interval.
For type '"prob"' a matrix with number of thresholds is produced. Each row corresponds to a row in newdata and contains the predicted probabilities to fall below the corresponding threshold.
For types '"class"', '"location"', and "'scale"' a vector is returned. 
For type '"class"' the interval with the highest probability is returned.
For types '"location"' and "'scale"' location and scale of the predicted latent distribution is returned respectively.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Jakob W. Messner
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hxlr}}
}
%\examples{

%}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~regression }
