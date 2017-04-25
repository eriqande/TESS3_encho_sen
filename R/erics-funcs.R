

#' This function returns the rasterStack that underlies the plots that you would get
#' if you were to send the same arguments to the function \code{\link{plot.tess3Q}},
#' 'cuz it might be nice to get at the underlying plot data.
#' @title Return rasters underlying the geographic maps of ancestry coefficients
#' @author Kevin Caye, Flora Jay, Olivier François (repackaged by Eric C. Anderson)
#'
#' @param coord a numeric matrix of dimension \eqn{n} times 2 where \eqn{n} is the
#' number of individuals. The matrix must contain (Longitude, Latitude) coordinates for all individuals.
#' @param resolution an integer vector \code{resolution = c(rx,ry)} for the resolution of the grid used to
#' computed the interpolating surface. \code{rx} and \code{ry} are resolution numbers for the x-axis and y-axis respectively.
#' @param window an integer vector for the plotting window, such that \code{window = c(xmin, xmax, ymin, ymax)}
#' contains the window's min and max coordinates.
#' @param background if \code{TRUE} a raster file is used as a stencil to render only raster pixel on earth.
#' @param map.polygon an object of class \code{sp::SpatialPolygonsDataFrame} used to crop the interpolating surface.
#' If \code{NULL}, the function \code{\link[rworldmap]{getMap}} is used.
#' @param raster.filename a raster file name used to compute the background stencil.
#' This is an alternative method to crop the interpolating surfaces. The default method uses \code{map.polygon}.
#' @param col.palette a list of color palettes. Color palettes can be defined by using the function \code{\link{CreatePalette}}.
#' @param interpolation.model an interpolation model used to compute the interpolating surface. Interpolation models can use the
#' functions \code{\link{FieldsTpsModel}} or \code{\link{FieldsKrigModel}}.
#' @param x an object of class \code{tess3Q} containing an ancestry coefficient matrix computed from \code{\link{tess3}}.
#' @param method THIS DOESN'T HAVE ANY EFFECT HERE...THE rasters that come back are the same.  a character string \code{"map.max"} or \code{"map.all"}. If \code{"map.all"}, interpolating surfaces are displayed
#' for the ancestry coefficients of all populations. If \code{"map.max"} the union of interpolating surfaces is displayed.
#' @param ... \code{\link{plot.default}} other parameters.
#'
#' @return None
#' @examples
#' library(tess3r)
#'
#' # Retrieve a dataset
#' data(data.at)
#'
#' # Run of TESS3
#' obj <- tess3(data.at$X, coord = data.at$coord, K = 5,
#'                  ploidy = 1, method = "projected.ls", openMP.core.num = 4)
#'
#' # Get the ancestry matrix
#' Q.matrix <- qmatrix(obj, K = 5)
#'
#' # Plot the spatial interpolation of the ancestry matrix
#' boing <- tess3Q_map_rasters(Q.matrix, data.at$coord, method = "map.max",
#'      resolution = c(400,400),
#'      interpolation.model = FieldsKrigModel(10), cex = .4,
#'      xlab = "Longitude", ylab= "Latitude", main = "Ancestry coefficients")
#' @export
tess3Q_map_rasters <- function(x, coord,
                        method = "map.max",
                        resolution = c(300,300), window = NULL,
                        background = TRUE, map.polygon = NULL,
                        raster.filename = NULL,
                        interpolation.model = FieldsKrigModel(10),
                        col.palette = CreatePalette(), ...) {

  CheckPlotParam(x, coord,
                 method,
                 resolution, window,
                 background, map.polygon,
                 raster.filename,
                 interpolation.model,
                 col.palette)


  if (method == "piechart") {
    stop("In development.")
    # PlotPiechartAncestryCoef(x, coord, window, background, col, ...)
  } else {
    interpol.stack <- ComputeInterpolStack(x, coord,
                                           interpolation.model, resolution, window,
                                           background, map.polygon, raster.filename)
  }
  return(interpol.stack)
}
