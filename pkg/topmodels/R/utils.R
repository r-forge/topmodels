#my_adjust_transparency <- function(col, alpha, default = 1) {
#  if (!isTRUE(alpha)) {
#    return(colorspace::adjust_transparency(col = col, alpha = alpha))
#  } else {
#    str_alpha <- substr(colorspace::adjust_transparency(col, alpha = TRUE), 8, 9)
#    str_col <- substr(colorspace::adjust_transparency(col, alpha = TRUE), 1, 7)
#    if (str_alpha == "FF") {
#      return(colorspace::adjust_transparency(str_col, alpha = default))
#    } else { 
#      return(col)
#    }
#  }
#}

set_minimum_transparency <- function(col, alpha_min) {
  ## sanity checks
  ## `col` w/i `extract_transparency()`
  stopifnot(
    is.numeric(alpha_min),
    length(alpha_min) == 1
  )

  ## get alpha 
  alpha <- extract_transparency(col, mode = "numeric", default = alpha_min)

  ## set alpha (minimum if necessary)
  alpha <- ifelse(alpha > alpha_min, alpha_min, alpha)

  col <- adjust_transparency(col, alpha = alpha)
 
  return(col) 
}

## FIXME: (ML) Use colorspace version once published on CRAN
adjust_transparency <- function(col, alpha = TRUE) {

  ## alpha argument controls new alpha values
  new_alpha <- alpha
  if(inherits(new_alpha, "hexmode")) new_alpha <- unclass(new_alpha)

  ## support S4 color specifications as well (with default alpha = 1)
  if(inherits(col, "color")) col <- paste0(colorspace::hex(col), "FF")

  ## keep indizes of NA colors
  NAidx <- which(is.na(col))
  n <- if(is.matrix(col) && is.numeric(col)) NCOL(col) else length(col)

  ## col has to be hex code, otherwise col2rgb is used
  if(is.character(col) &&
    (all(substr(col, 1L, 1L) == "#") & all(nchar(col) %in% c(7L, 9L))))
  {
    ## extract alpha from hex (if any)
    alpha <- substr(col, 8L, 9L)
    ## retain only RGB in hex
    col <- substr(col, 1L, 7L)
  } else {
    if(!(is.matrix(col) && is.numeric(col))) col <- col2rgb(col, alpha = TRUE)
    ## extract alpha values (if any)
    alpha <- if(NROW(col) > 3L) convert_transparency(col[4L, ]/255, mode = "character") else rep.int("FF", n)
    ## retain only RGB
    col <- rgb(col[1L, ], col[2L, ], col[3L, ], maxColorValue = 255)
  }

  ## adjust alpha transparency
  if(is.null(new_alpha)) {
    alpha[alpha == "FF"] <- ""
  } else if(identical(new_alpha, FALSE)) {
    alpha <- rep.int("", n)
  } else if(identical(new_alpha, TRUE)) {
    alpha[alpha == ""] <- "FF"
  } else {
    alpha <- convert_transparency(new_alpha, mode = "character")
    if(length(alpha) != n) {
      n <- max(n, length(alpha))
      col <- rep_len(col, n)
      alpha <- rep_len(alpha, n)
    }
  }

  ## add alpha again (if any) and manage NAs
  col <- paste(col, alpha, sep = "")
  if(length(NAidx) > 0) col[NAidx] <- NA

  return(col)
}

## FIXME: (ML) Use colorspace version once published on CRAN
extract_transparency <- function(col, mode = "numeric", default = 1) {

  ## handle mode of return value
  mode <- match.arg(mode, c("numeric", "double", "integer", "character", "hexmode"))
  if(mode == "double") mode <- "numeric"

  ## handle default
  if(length(default) > 1L) {
    warning("'default' should be of length 1, first element used")
    default <- default[1L]
  } else if(length(default) < 1L) {
    default <- NA
  }
  if(inherits(default, "hexmode")) default <- unclass(default)
  if(mode == "hexmode") {
    hexmode <- TRUE
    mode <- "integer"
  } else {
    hexmode <- FALSE
  }
  na <- switch(mode,
    "numeric" = NA_real_,
    "character" = NA_character_,
    NA_integer_)
  if(is.na(default)) {
    default <- na
  } else if(is.character(default) | is.numeric(default)) {
    default <- convert_transparency(default, mode = mode)
  } else {
    warning("unknown type of 'default' using NA instead")
    default <- na
  }

  ## for S4 color specifications or RGB matrices the default is used
  if(inherits(col, "color")) col <- colorspace::coords(col)

  ## number of colors and position of NA colors
  if(is.matrix(col) && is.numeric(col)) {
    n <- NCOL(col)
    ina <- which(apply(is.na(col), 1L, any))
  } else {
    n <- length(col)
    ina <- which(is.na(col))
  }

  ## cases where alpha can be extracted
  ialpha <- if(is.character(col)) {
    which(substr(col, 1L, 1L) == "#" & nchar(col) == 9L)
  } else {
    integer(0)
  }

  ## set up return value
  alpha <- rep.int(default, n)
  if(length(ialpha) > 0L) alpha[ialpha] <- convert_transparency(substr(col[ialpha], 8L, 9L), mode = mode)
  alpha[ina] <- na
  if(hexmode) alpha <- structure(alpha, class = "hexmode")
  return(alpha)
}

## FIXME: (ML) Use colorspace version once published on CRAN
convert_transparency <- function(x, mode = "numeric") {
  mode <- match.arg(mode, c("numeric", "integer", "character"))
  if(is.character(x)) {
    if(!all(toupper(x) %in% format(as.hexmode(0L:255L), width = 2L, upper.case = TRUE))) {
      stop("invalid character specification of alpha transparency, must be in 00, 01, ..., FF")
    }
    x <- switch(mode,
      "integer" = as.integer(as.hexmode(x)),
      "numeric" = as.integer(as.hexmode(x))/255,
      x)
  } else if(is.integer(x)) {
    if(!all(x %in% 0L:255L)) {
      stop("invalid integer specification of alpha transparency, must be in 0L, 1L, ..., 255L")
    }
    x <- switch(mode,
      "numeric" = x/255,
      "character" = format(as.hexmode(x), width = 2L, upper.case = TRUE),
      x)
  } else if(is.numeric(x)) {
    if(!all(x >= 0 & x <= 1)) {
      stop("invalid numeric specification of alpha transparency, must be in [0, 1]")
    }
    x <- switch(mode,
      "integer" = as.integer(round(x * 255 + 0.0001)),
      "character" = format(as.hexmode(round(x * 255 + 0.0001)), width = 2L, upper.case = TRUE),
      x)
  }
  return(x)
}


annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
  {
    ggplot2::layer(data = data, stat = ggplot2::StatIdentity, position = ggplot2::PositionIdentity, 
          geom = ggplot2:::GeomCustomAnn,
          inherit.aes = TRUE, params = list(grob = grob, 
                                            xmin = xmin, xmax = xmax, 
                                            ymin = ymin, ymax = ymax))
  }
