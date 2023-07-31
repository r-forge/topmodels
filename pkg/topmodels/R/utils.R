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

#' Helper functions to check available support for some S3methods
#'
#' @param methods character, name of the method (e.g., is_continuous, print, ...)
#' @classes character vector of length > 0, classes to check.
#'
#' @return Returns TRUE if the method exists for one of the given classes, else FALSE.
#'
#' @importFrom utils getS3method
hasS3method <- function(method, classes) {
  any(sapply(classes, function(cls) {
    tryCatch(is.function(getS3method(method, class = cls)), error = function(e) FALSE)
  }))
}

annotation_custom2 <- function(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) {
  ggplot2::layer(
    data = data, stat = ggplot2::StatIdentity, position = ggplot2::PositionIdentity,
    geom = ggplot2::GeomCustomAnn,
    inherit.aes = TRUE, params = list(
      grob = grob,
      xmin = xmin, xmax = xmax,
      ymin = ymin, ymax = ymax
    )
  )
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
  stopifnot(is.null(default) || length(default) >= 1)

  ## helper_function
  is_fun_or_unique <- function(x) {
    if (is.function(x) || length(unique(x)) == 1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }

  ## use arg if defined by user
  arg_by_user <- try(eval(parse(text = sprintf("!missing(%s)", arg_name)), envir = parent.frame()), 
    silent = TRUE)
  if (!inherits(arg_by_user, "try-error") && isTRUE(arg_by_user)) {
    arg_fun <- eval(parse(text = arg_name), envir)
  } else {
    arg_fun <- NULL
  }

  ## get arg value from attributes
  arg_attr <- attr(object, arg_name)

  ## conditional return
  if (force_single) {
    if (is.null(arg_fun) && is.null(arg_attr)) {
      rval <- default[1L]
    } else if (is.null(arg_fun) && !is_fun_or_unique(arg_attr)) {
      message(sprintf(
        " * as arg `%s`'s definition is not unique w/i object's attributes, using the default (\"%s\")",
        arg_name, default[1L]
      ))
      rval <- default[1L]
    } else if (is.null(arg_fun)) {
      rval <- if (is.function(arg_attr)) arg_attr else arg_attr[[1]]
    } else {
      rval <- if (is.function(arg_fun)) arg_fun else arg_fun[[1]]
    }
  } else {
    if (is.null(arg_fun) && is.null(arg_attr)) {
      rval <- default
    } else if (is.null(arg_fun)) {
      rval <- arg_attr
    } else {
      rval <- arg_fun
    }
  }
  
  return(rval)
}


prepare_arg_for_attributes <- function(object,
                                       arg_name,
                                       missing = NA,
                                       force_single = FALSE) {
  rval <- unlist(
    lapply(
      object,
      function(x) ifelse(is.null(attr(x, arg_name)), missing, attr(x, arg_name))
    )
  )

  if (force_single) {
    if (length(unique(rval)) > 1) {
      message(sprintf(
        " * as arg `%s`'s definition is not unique, using solely the first (\"%s\")", arg_name, rval[1]
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


compute_breaks <- function(mid, width, offset = FALSE) {

  ## sanity checks
  stopifnot(is.numeric(mid), is.null(dim(mid)))
  stopifnot(is.numeric(width), is.null(dim(width)))
  stopifnot(length(mid) == length(width))
  stopifnot(is.logical(offset))

  ## compute breaks
  if (!offset) { 
    n <- length(width)
    rval <- c(mid - width / 2, mid[n] + width[n] / 2)
  } else {
    xleft <- mid - width / 2
    xright <- mid + width / 2
    rval <- c( 
      xleft[1],
      colMeans(rbind(xleft[-1], xright[-length(xright)])),
      xright[length(xright)]
    )
  }
  rval
}

set_aes_helper_geoms <- function(aes_geom, aes_arg, aes_default = NULL) {

  stopifnot(is.list(aes_geom), is.list(aes_arg), is.null(aes_default) || is.list(aes_default))

  for (i in names(aes_geom)) {
    if (is.null(aes_arg[[i]]) && !is.null(aes_default[[i]])) {
      aes_geom[[i]] <- aes_default[[i]]
    } else if (!is.null(aes_arg[[i]])) {
      aes_geom[[i]] <- aes_arg[[i]]
    } 
  }
  aes_geom
}

# Convert Wide to Long Format
#
# Converting wide \code{data.frame} to long format for the functionality
# needed in topmodels (base R).
#
# @param x data.frame
# @param id_cols character vector of length >= 1; columns to be kept as is, used
#        for grouping.
# @param keep_cols \code{NULL} or character, columns to keep (as is).
# @param values_from character of length 1, name of the column from which the
#        values should be taken.
# @param names_to character of length 1, name of the variable to store the original
#        variable name (\code{values_from}).
# @param values_to character of length 1, name of the variable to store the data.
# @param check logical, if \code{TRUE} a series of sanity checks are performed.
#
# @return Returns the reshaped \code{data.frame}.
#
# @author Reto
wide_to_long <- function(x, id_cols, keep_cols, values_from, names_to, values_to, check = FALSE) {
    stopifnot(isTRUE(check) || isFALSE(check))
    if (check) {
        stopifnot(is.data.frame(x))
        stopifnot(is.character(values_from))
        stopifnot(is.character(names_to) && length(names_to) == 1)
        stopifnot(is.character(values_to) && length(values_to) == 1)
        stopifnot(all(id_cols %in% names(x)), values_from %in% names(x))
        stopifnot(is.null(keep_cols) || is.character(keep_cols))
        if (is.character(keep_cols)) stopifnot(all(keep_cols %in% names(x)))
    }

    # Splitting the data
    x <- split(x[c(id_cols, keep_cols, values_from)], interaction(x[id_cols], drop = TRUE))
    # Scoping various elements
    fn <- function(y, vf) {
        y <- cbind(y[, c(id_cols, keep_cols), drop = FALSE], vf, y[, vf])
        names(y)[ncol(y) - 1:0] <- c(names_to, values_to)
        return(y)
    }
    x <- do.call(rbind, lapply(values_from, function(vf) do.call(rbind, lapply(x, fn, vf = vf))))
    return(structure(transform(x, group = as.character(group)),
                     row.names = seq_len(nrow(x)))) # Reset row names
}
