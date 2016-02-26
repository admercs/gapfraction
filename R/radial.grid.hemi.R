#' Modified Radial Grid Function for Hemispherical Lens Geometries
#'
#' This function creates radial grid plots for four common hemispherical lens geometries: equi-distant, equi-angular, stereographic, and orthographic
#' @param labels Parameter adopted from the radial.grid function of the plotrix package. Defaults to NA.
#' @param label.pos Parameter adopted from the radial.grid function of the plotrix package. Defaults to NULL.
#' @param radlab Parameter adopted from the radial.grid function of the plotrix package. Defaults to FALSE.
#' @param radial.lim Parameter adopted from the radial.grid function of the plotrix package. Defaults to NULL.
#' @param start Parameter adopted from the radial.grid function of the plotrix package. Defaults to 0.
#' @param clockwise Parameter adopted from the radial.grid function of the plotrix package. Defaults to FALSE.
#' @param label.prop Parameter adopted from the radial.grid function of the plotrix package. Defaults to 1.1.
#' @param grid.pos Parameter adopted from the radial.grid function of the plotrix package. Defaults to grid.pos.
#' @param grid.col Parameter adopted from the radial.grid function of the plotrix package. Defaults to gray.
#' @param grid.bg Parameter adopted from the radial.grid function of the plotrix package. Defaults to transparent.
#' @param show.radial.grid Parameter adopted from the radial.grid function of the plotrix package. Defaults to TRUE.
#' @param model Hemispherical lens geometry model to use. Options include equi-distant (\code{equidist}), equi-angular (\code{equiangle}), stereographic (\code{stereo}), and orthographic (\code{ortho}). Defaults to NA.
#' @param r Hemispherical lens geometry string to use if no \code{model} is specified. String functions for \code{r} must be of the form \code{r = "tan(theta/2)"}. Defaults to theta.
#' @keywords radial, grid, hemispherical
#' @export
#' @return The results of \code{radial.grid.hemi}
#' @examples
#' radial.grid.hemi()

radial.grid.hemi <- function (labels = NA, label.pos = NULL, radlab = FALSE, radial.lim = NULL, start = 0,
                              clockwise = FALSE, label.prop = 1.1, grid.pos = grid.pos, grid.col = "gray",
                              grid.bg = "transparent", show.radial.grid = TRUE, model = NA, r = "theta")
{
  if(model == "stereo")    r <- "tan(theta/2)"
  if(model == "ortho")     r <- "sin(theta)"
  if(model == "equidist")  r <- "theta"
  if(model == "equiangle") r <- "sin(theta/2)"

  rho   <- gsub("theta", "grid.pos[i]", r)
  r.phi <- gsub("theta", "(pi/2)", r)

  par(xpd = TRUE)
  if (is.null(label.pos))
    label.pos <- seq(0, 1.8 * pi, length = 9)
  if (!is.null(labels)) {
    if (is.na(labels[1]))
      labels <- as.character(round(label.pos, 2))
  }
  if (clockwise)
    label.pos <- -label.pos
  if (start)
    label.pos <- label.pos + start
  phi <- seq(0, 1.96 * pi, by = 0.04 * pi)
  for (i in seq(length(grid.pos), 1, by = -1)) {
    xpos <- cos(phi) * (eval(parse(text=rho)) - radial.lim[1])
    ypos <- sin(phi) * (eval(parse(text=rho)) - radial.lim[1])
    polygon(xpos, ypos, border = grid.col, col = grid.bg)
  }
  maxlength <- max(grid.pos) - radial.lim[1]
  if (show.radial.grid) {
    xpos <- cos(label.pos) * eval(parse(text=r.phi)) #maxlength
    ypos <- sin(label.pos) * eval(parse(text=r.phi)) #maxlength
    segments(0, 0, xpos, ypos, col = grid.col)
    xpos <- cos(label.pos) * eval(parse(text=r.phi)) #maxlength
    ypos <- sin(label.pos) * eval(parse(text=r.phi)) #maxlength
  }
  if (!is.null(labels)) {
    xpos <- cos(label.pos) * label.prop * eval(parse(text=r.phi)) #maxlength
    ypos <- sin(label.pos) * label.prop * eval(parse(text=r.phi)) #maxlength
    if (radlab) {
      for (label in 1:length(labels)) {
        if (radlab < 0)
          labelsrt <- 180 * label.pos[label]/pi - 90 + 180 * (label.pos[label] > pi && label.pos[label] < 2 * pi)
        else labelsrt <- (180 * label.pos[label]/pi) + 180 * (label.pos[label] > pi/2 && label.pos[label] < 3 * pi/2)
        text(xpos[label], ypos[label], labels[label], cex = par("cex.axis"), srt = labelsrt)
      }
    }
    else plotrix::boxed.labels(xpos, ypos, labels, ypad = 0.7, border = FALSE, cex = par("cex.axis"))
  }
  par(xpd = FALSE)
}
