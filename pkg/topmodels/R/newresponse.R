newresponse <- function(object, ...) {
  UseMethod("newresponse")
}

newresponse.default <- function(object, newdata, na.action = na.pass, ...)
{
  if(missing(newdata) || is.null(newdata)) {
    model.response(model.frame(object))
  } else {
    model.response(model.frame(update(terms(object), . ~ 1), newdata, na.action = na.action, ...))
  }
}
