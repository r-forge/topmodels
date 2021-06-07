topmodels <- function(object,
                      plot = TRUE,
                      class = NULL,
                      newdata = NULL,
                      na.action = na.pass,
                      which = NULL,
                      ask = dev.interactive(), # FIXME: (ML) Does not work for ggplot.
                      spar = TRUE, # FIXME: (ML) What does this do? Needed? Does not work for ggplot.
                      single_page = NULL, # FIXME: (ML) Does not work for ggplot.
                      ...) {
  # -------------------------------------------------------------------
  # SET UP PRELIMINARIES
  # -------------------------------------------------------------------
  ## guess plotting flavor
  if (isFALSE(plot)) {
    plot <- "none"
  } else if (isTRUE(plot)) {
    plot <- if ("package:ggplot2" %in% search()) "ggplot2" else "base"
  } else if (!is.character(plot)) {
    plot <- "base"
  }
  plot <- try(match.arg(plot, c("none", "base", "ggplot2")))
  stopifnot(
    "`plot` must either be logical, or match the arguments 'none', 'base' or 'ggplot2'" =
      !inherits(plot, "try-error")
  )

  ## guess output class
  if (is.null(class)) {
    class <- if ("package:tibble" %in% search()) "tibble" else "data.frame"
  }
  class <- try(match.arg(class, c("tibble", "data.frame")))
  stopifnot(
    "`class` must either be NULL, or match the arguments 'tibble', or 'data.frame'" =
      !inherits(class, "try-error")
  )

  ## check if S3 methods exist
  if (!any(class(object) %in% gsub("procast.", "", methods("procast")))) {
    stop(
      sprintf(
        "The `object` must be one of the following classes: %s.",
        paste(gsub("procast.", "", methods("procast")), collapse = ", ")
      )
    )
  }

  # check which method should be plotted
  which.match <- c("rootogram", "pithist", "reliagram", "qqrplot", "wormplot")
  if (is.null(which)) which <- which.match
  if (!is.character(which)) {
    if (any(which > 5L)) {
      which <- which[which <= 5L]
    }
    which <- which.match[which]
  } else {
    which <- which.match[pmatch(tolower(which), which.match)]
  }

  if (length(which) > length(which.match) || !any(which %in% which.match)) {
    stop("argument which is specified wrong!")
  }

  ## define layout of plot
  if (spar && !plot == "ggplot2") {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
  }

  if (isTRUE(single_page)) {
    ask <- FALSE
  }

  if (prod(par("mfcol")) >= length(which)) {
    spar <- FALSE
    ask <- FALSE
  }

  if (spar && !plot == "ggplot2") {
    if (!ask) {
      par(mfrow = n2mfrow(length(which)))
    } else {
      par(ask = ask)
    }
  }

  # -------------------------------------------------------------------
  # CHECK ARGUMENTS AND PREPARE FUNCTION CALLS
  # -------------------------------------------------------------------
  ## get call and set topmodel() args to NULL
  mc <- as.list(match.call())[-1]
  mc$which <- NULL
  mc$ask <- NULL
  mc$spar <- NULL
  mc$single_page <- NULL

  ## get possible arguments
  arg_avail <- list()

  if (plot == "base") {
    arg_avail[["rootogram"]] <- unique(names(c(formals(rootogram.default), formals(plot.rootogram))))
    arg_avail[["pithist"]] <- unique(names(c(formals(pithist.default), formals(plot.pithist))))
    arg_avail[["reliagram"]] <- unique(names(c(formals(reliagram.default), formals(plot.reliagram))))
    arg_avail[["qqrplot"]] <- unique(names(c(formals(qqrplot.default), formals(plot.qqrplot))))
    arg_avail[["wormplot"]] <- unique(names(c(formals(wormplot.default), formals(plot.wormplot))))
  } else {
    arg_avail[["rootogram"]] <- unique(names(c(formals(rootogram.default), formals(autoplot.rootogram))))
    arg_avail[["pithist"]] <- unique(names(c(formals(pithist.default), formals(autoplot.pithist))))
    arg_avail[["reliagram"]] <- unique(names(c(formals(reliagram.default), formals(autoplot.reliagram))))
    arg_avail[["qqrplot"]] <- unique(names(c(formals(qqrplot.default), formals(autoplot.qqrplot))))
    arg_avail[["wormplot"]] <- unique(names(c(formals(wormplot.default), formals(autoplot.wormplot))))
  }

  arg_avail <- lapply(arg_avail, function(x) unique(c(x, names(par()))))
  arg_avail <- lapply(arg_avail, function(x) x[x != "..."])

  ## warning if any provided arg is not valid argument
  if (any(!(names(mc) %in% unique(do.call("c", arg_avail))))) {
    warning(sprintf(
      "The following arguments are not valid in any of the topmodels' plotting functions: %s",
      paste0(names(mc)[!(names(mc) %in% unique(do.call("c", arg_avail)))], collapse = ",")
    ))
  }

  ## prepare a list of arguments for each plot and fill it with arguments
  mc <- lapply(mc, function(x) if (typeof(x) == "language") eval(x) else x)

  ## remove all invalid names
  tmp_check <- sapply(mc, function(x) (length(x) > 1 && is.null(names(x))) || any(!names(x) %in% which))
  if (any(tmp_check)) {
    warning(sprintf(
      "The following arguments are not correctly specified and therefore (partially) omitted: %s",
      paste0(names(mc)[tmp_check], collapse = ",")
    ))
  }

  ## check if named vector is provided and any of the names matches function name
  arg <- list("rootogram" = mc, "pithist" = mc, "reliagram" = mc, "qqrplot" = mc, "wormplot" = mc)
  arg <- lapply(
    seq_along(arg),
    function(idx) {
      lapply(
        arg[[idx]],
        function(x) {
          if (is.null(names(x))) {
            x
          } else if (!is.null(names(x)) && !names(arg)[idx] %in% names(x)) {
            NULL
          } else {
            x[names(x) == names(arg)[idx]]
          }
        }
      )
    }
  )

  ## remove NULLs from list
  arg <- lapply(arg, function(x) x[!sapply(x, is.null)])

  ## add proper titles
  arg <- lapply(seq_along(arg), function(idx) {
    arg[[idx]]$main <- c(
      "Rootogram", "PIT histogram", "Reliability diagram",
      "Q-Q residuals plot", "Worm plot"
    )[idx]
    arg[[idx]]
  })

  ## look up if provided arguments match possible arguments
  arg <- lapply(seq_along(arg), function(idx) arg[[idx]][names(arg[[idx]]) %in% arg_avail[[idx]]])

  ## set names
  names(arg) <- c("rootogram", "pithist", "reliagram", "qqrplot", "wormplot")


  # -------------------------------------------------------------------
  # FUNCTION CALLS: CALCULATE AND PLOT
  # -------------------------------------------------------------------
  rval <- list()
  if (plot == "base") {
    ## calculate and plot
    if ("rootogram" %in% which) rval$rootogram <- do.call(rootogram, arg[[1]])
    if ("pithist" %in% which) rval$pithist <- do.call(pithist, arg[[2]])
    if ("reliagram" %in% which) rval$reliagram <- do.call(reliagram, arg[[3]])
    if ("qqrplot" %in% which) rval$qqrplot <- do.call(qqrplot, arg[[4]])
    if ("wormplot" %in% which) rval$wormplot <- do.call(wormplot, arg[[5]])
  } else {
    ## first calculate
    arg <- lapply(arg, function(x) {
      x$plot <- FALSE
      x
    })
    if ("rootogram" %in% which) rval$rootogram <- do.call(rootogram, arg[[1]])
    if ("pithist" %in% which) rval$pithist <- do.call(pithist, arg[[2]])
    if ("reliagram" %in% which) rval$reliagram <- do.call(reliagram, arg[[3]])
    if ("qqrplot" %in% which) rval$qqrplot <- do.call(qqrplot, arg[[4]])
    if ("wormplot" %in% which) rval$wormplot <- do.call(wormplot, arg[[5]])

    ## use return val as new object
    for (name in names(rval)) {
      arg[[name]]$object <- rval[[name]]
    }

    ## get plotting dimension for grid
    ## TODO: (ML) Can you specify layout.pos.col w/o exact location? Can maybe be improved.
    plt_dim <- n2mfrow(length(which))
    plt_names <- names(rval)
    plt_rows <- rep_len(rep(1:plt_dim[1], each = plt_dim[2]), length(plt_names))
    plt_cols <- rep_len(1:plt_dim[2], length(plt_names))

    ## set up grid and plot
    ## TODO: (ML)
    ## * Works as additional args of main funs are not used in autoplot (except `names(par())`)
    ## * Rewrite cleaner version.
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(plt_dim[1], plt_dim[2])))
    for (idx in seq_along(plt_names)) {
      print(
        do.call(ggplot2::autoplot, arg[[plt_names[idx]]]),
        vp = grid::viewport(layout.pos.row = plt_rows[idx], layout.pos.col = plt_cols[idx])
      )
    }
  }

  # -------------------------------------------------------------------
  # RETURN INVISIBLY
  # -------------------------------------------------------------------
  return(invisible(rval))
}
