## TODO:
## - Z: create R-Forge repos "topmodels" [done]
## - Z: svndump "crch" -> import to "topmodels" on R-Forge [done]
## - Z: package skeleton "topmodels" on R-Forge [done]
## - Z: procast() generic plus flexible procast_setup() [done]
## - JM: procast.crch() with new procast_setup() infrastructure.
## - JM: reliagram() prototype
## - JM: port pithist() from "countreg" + add type = "proportional" + Manu style
## - Z: port qqrplot() and rootogram() from "countreg" plus adaptations
## - JM/Z/.../CK/IK/.../SW/MS/GS/???: procast methods for betareg, countreg, glm, gamlss, mgcv::gam, ...

## TODO:
## - scoring rules for fitted model objects:
##   in-sample vs. out-of-sample / aggregated vs. observation-wise contributions
## - out-of-sample logLik()/logs(), crps(), ..., discretized log-score


topmodels <- function(object, 
                      newdata = NULL,
                      na.action = na.pass, 
                      which = NULL,
                      ask = dev.interactive(),
                      spar = TRUE,
                      pages = NULL,
                      ...) {

  ## check if S3 methods exist
  if (!any(class(object) %in% gsub("procast.", "", methods("procast")))) {
    stop(
      sprintf("The `object` must be one of the following classes: %s.", 
      paste(gsub("procast.", "", methods("procast")), collapse = ", "))
    )
  }

  # check which method should be plotted
  which.match <- c("pithist", "qqrplot", "reliagram", "rootogram", "wormplot")
  if (is.null(which)) which <- which.match
  if (!is.character(which)) {
    if (any(which > 5L))
      which <- which[which <= 5L]
    which <- which.match[which]
  } else {
    which <- which.match[pmatch(tolower(which), which.match)]
  }

  if (length(which) > length(which.match) || !any(which %in% which.match)) {
    stop("argument which is specified wrong!")
  }


  ## define layout of plot
  if (spar) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
  }

  if (!is.null(pages)) {
    ask <- !(pages == 1)
  }

  if (prod(par("mfcol")) > 1L) {
    spar <- FALSE
    ask <- FALSE
  }

  if (spar) {
    if (!ask) {
      par(mfrow = n2mfrow(length(which)))
    } else {
      par(ask = ask)
    }
  }

  ## prepare plotting arguments (needed to display main titles correct)
  mc <- as.list(match.call())[-1]
  mc$which <- NULL
  mc$ask <- NULL
  mc$spar <- NULL
  mc$pages <- NULL
  # FIXME: (ML) Make only valid arguments in ... are passed on to functions

  rval <- list()
  if ("rootogram" %in% which) rval$rootogram <- do.call(rootogram, mc)
  if ("pithist"   %in% which) rval$pithist <-   do.call(pithist, mc)
  if ("reliagram" %in% which) rval$reliagram <- do.call(reliagram, mc)
  if ("qqrplot"   %in% which) rval$qqrplot <-   do.call(qqrplot, mc)
  if ("wormplot"  %in% which) rval$wormplot <-  do.call(wormplot, mc)

  return(invisible(rval))
}

                       

  
