set_minimum_transparency <- function(col, alpha_min) {
  ## sanity checks
  ## `col` w/i `extract_transparency()`
  stopifnot(
    is.numeric(alpha_min),
    length(alpha_min) == 1
  )

  ## get alpha 
  alpha <- colorspace::extract_transparency(col, mode = "numeric", default = alpha_min)

  ## set alpha (minimum if necessary)
  alpha <- ifelse(alpha > alpha_min, alpha_min, alpha)

  col <- colorspace::adjust_transparency(col, alpha = alpha)
 
  return(col) 
}


annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) {
    ggplot2::layer(data = data, stat = ggplot2::StatIdentity, position = ggplot2::PositionIdentity, 
          geom = ggplot2::GeomCustomAnn,
          inherit.aes = TRUE, params = list(grob = grob, 
                                            xmin = xmin, xmax = xmax, 
                                            ymin = ymin, ymax = ymax))
}


# Helper function inspired by internal from `ggplot2` defined in `performance.R`
# (needed for geom_*s to modify aesthetics for different styles)
my_modify_list <- function(old, new, force = FALSE) {

  if (force) {
    for (i in names(new)) old[[i]] <- new[[i]]
  } else {
    for (i in names(new)) old[[i]] <- if (all(is.na(old[[i]]) | old[[i]] == 0.999)) new[[i]] else old[[i]]
  }

  old
}


## Helper function for setting arguments to values w/i attributes
use_arg_from_attributes <- function(object, 
                                    arg_name, 
                                    default = NULL, 
                                    force_single = FALSE, 
                                    envir = parent.frame()) {

  ## check input     
  stopifnot(is.character(arg_name), length(arg_name) == 1)

  ## get arg value from function
  arg_fun <- try(eval(parse(text = arg_name), envir))
  if (inherits(arg_fun, "try-error")) stop(sprintf("arg `%s` is not defined", arg_name))

  ## get arg value from attributes
  arg_attr <- attr(object, arg_name)

  ## conditional return 
  if (is.null(arg_fun) && force_single && length(unique(arg_attr)) > 1) {
    message(sprintf(
      " * as arg `%s`'s definition is not unique w/i object's attributes, using the default", arg_name
    ))
    rval <- default
  } else if (is.null(arg_fun) && is.null(arg_attr)) {
    rval <- default
  } else if (is.null(arg_fun)){
    rval <- arg_attr
  } else {
    rval <- arg_fun
  }

  if (force_single) {
    if (length(unique(rval)) > 1) {
      message(sprintf(
        " * as arg `%s`'s definition is not unique, using solely the first", arg_name
      ))
    }
    return(rval[1])
  } else {
    return(rval)
  }
}


duplicate_last_value <- function(x) {
  
  ## sanity checks
  stopifnot(is.null(dim(x)))
 
  ## repeat last value and return  
  c(x, x[NROW(x)])
}

