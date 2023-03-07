# Code to do image correlation to help adjust plot locations
#

# isOpt
#
#'
#' This is a helper function used in the package to test the
#' value of an object: usually a function parameter. While this simply
#' wraps \code{is.null()}, it makes for slightly more readable code when
#' checking function parameters.
#'
#' @param x object
#' @return A (invisible) boolean value, TRUE indicates the object value is not NULL,
#'   FALSE indicates the object value is NULL.
#' @keywords internal
isOpt <- function(x) { invisible(!is.null(x)) }

#' computeCorrelationAndHeightError
#'
#' Computes correlation between the stem map (with height) and CHM along with a
#' height error statistic.
#'
#' This function is typically used in a brute force process designed to find the
#' best location for a plot given stem-mapped trees and a lidar-derived canopy
#' height model (CHM). It can also be used to compute a one-time set of fit
#' metrics useful for comparing one location with another.
#'
#' The correlation is Pearson's correlation coefficient computed by
#' terra::layerCor between a tree map and the CHM (masked by the tree map raster
#' layer). The tree map is created using a buffer around the tree locations
#' (buffer can be 0) with cells for the tree populated with height values. Using
#' a small buffer seems to help in situations where trees lean (base XY and top
#' XY are different) or where the top of the tree is not best represented by a
#' single pixel in the CHM (possibly deciduous trees). In classic applications
#' that use the correlation coefficient to compare images, gray scale values are
#' used to compute the coefficient. In our case, we use the height values from
#' the CHM and tree map to compute the coefficient.
#'
#' The height error is the weighted mean absolute error between the tree map and
#' the CHM.
#'
#' The is method is taken from:
#' Monnet, Jean-Matthieu, and Éric Mermin. 2014. "Cross-Correlation of Diameter Measures for
#' the Co-Registration of Forest Inventory Plots with Airborne Laser Scanning Data" Forests 5,
#' no. 9: 2307-2326. https://doi.org/10.3390/f5092307
#'
#' The application of Pearson's correlation coefficient is also described in:
#' Arthur Miranda Neto, Alessandro Corrêa Victorino, Isabelle Fantoni, Douglas Eduardo
#' Zampieri, Janito Vaqueiro Ferreira, et al.. Image Processing Using Pearson’s
#' Correlation Coefficient: Applications on Autonomous Robotics. 13th International
#' Conference on Mobile Robots and Competitions (Robotica 2013), Apr 2013, Lisbon,
#' Portugal. pp.14-19. ffhal-00860912f
#'
#' @param stemMap Data frame with stem map data for a single plot. Needs to have XY coordinates for
#'  trees and a height for each tree.
#' @param coords Names of fields in \code{stemMap} containing X and Y values for tree locations.
#' @param htField Name of field in \code{stemMap} containing the tree height. Units should be the same
#'  as the height units in the CHM.
#' @param plotX Plot easting.
#' @param plotY plot northing.
#' @param plotRadius Radius of the plot used when measuring trees.
#' @param crs EPSG code for the coordinate reference system.
#' @param treeBufferSize Size of the buffer added to each tree. Typically small but this seems to help when there
#'  is error in relative tree locations. \code{treeBufferSize} can be a list of values so you can use a different
#'  buffer size for each tree. This may be useful if you want to scale the area considered
#'  when computing correlation and height error to the tree height. if \code{treeBufferSize} is 0,
#'  single cells in the tree map are populated with the tree height (tree locations are considered as
#'  points).
#' @param CHM Canopy height model in terra SpatRaster format.
#' @param searchRadius Maximum distance to offset plot location. This is needed to crop the CHM to the full
#'  area that will be tested for a good location to ensure that correlation and height error metrics match
#'  those returned by \code{testPlotLocations}.
#'
#' @return Data frame with the correlation and height error statistics. The combined error is not meaningful
#'  for a single location since the height error term in the combined score would go to 0.0.
#' @export
#'
#' @examples
computeCorrelationAndHeightError <- function(
  stemMap,
  coords = NULL,
  htField = NULL,
  plotX = 0,
  plotY = 0,
  plotRadius = 17.68,
  crs = NULL,
  treeBufferSize = 1.0,
  CHM = NULL,
  searchRadius = 0
) {
  # check parameters
  if (!isOpt(coords) || !isOpt(htField) || !isOpt(crs) || !isOpt(CHM))
    stop("Invalid parameters")

  # sort trees from shortest to tallest. This forces the highest height values into the
  # PM raster when trees overlap.
  stemMap <- stemMap[order(stemMap[, htField]),]

  # crop CHM to plot area expanded by the maximum search radius
  plotCHM <- terra::crop(CHM
                  , terra::ext(plotX - plotRadius - searchRadius,
                               plotX + plotRadius + searchRadius,
                               plotY - plotRadius - searchRadius,
                               plotY + plotRadius + searchRadius)
                  , snap = "near")

  # create stem map object
  trees <- sf::st_as_sf(stemMap,
                    coords = coords,
                    remove = FALSE,
                    crs = crs)

  # buffer trees
  if (length(treeBufferSize) > 1) {
    treesBuf <- sf::st_buffer(trees, treeBufferSize)
  }
  else {
    if (treeBufferSize > 0.0) {
      treesBuf <- sf::st_buffer(trees, treeBufferSize)
    }
    else {
      treesbuf <- trees
    }
  }

  r <- terra::rast(terra::vect(treesBuf)
            , extent = terra::ext(plotCHM)
            , resolution = terra::res(plotCHM)[1]
  )

  PM <- terra::rasterize(terra::vect(treesBuf)
                  , r
                  , field = "T3Ht"
                  , extent = terra::ext(plotCHM)
                  , resolution = terra::res(plotCHM)[1]
  )

  # stack layers so we can compute correlation
  # masking removes all CHM area not associated with plot trees (lidar visible)
  stack <- c(PM, terra::mask(plotCHM, PM))

  # Compute correlation and height error
  correlation <- terra::layerCor(stack, fun = "pearson", na.rm = TRUE)
  correlation <- correlation$pearson[1,2]

  # compute height error
  htError <- terra::global((abs(PM - plotCHM) * PM * PM), "sum", na.rm = TRUE) / terra::global(PM * PM, "sum", na.rm = TRUE)
  htError <- htError$sum

  return(data.frame(correlation, htError))
}

#' testPlotLocations
#'
#' Apply a brute-force approach to test a grid of potential plot locations and
#' report the results in a data frame including the offsets from the original
#' plot location, correlation scores and overall height error associated with
#' the trial plot location. The results are used with findBestPlotLocation() to
#' find the best location based on the correlation or the height error results
#' or a combination of the two. For the combined rule, the height error is
#' normalized relative to the maximum height error and subtracted from 1.0. The
#' combined "score" is then the correlation * normalized height and the best
#' location is the one with the highest combined score.
#'
#' The grid of potential plot locations is created using the resolution of the CHM so the
#' plotRadius should be a multiple of the CHM resolution.
#'
#' @param stemMap Data frame with stem map data for a single plot. Needs to have XY coordinates for
#'  trees and a height for each tree.
#' @param coords Names of fields in stemMap containing X and Y values for tree locations.
#' @param htField Name of field in stemMap containing the tree height. Units should be the same
#'  as the height units in the CHM.
#' @param plotX Plot easting.
#' @param plotY plot northing.
#' @param plotRadius Radius of the plot used when measuring trees.
#' @param crs EPSG code for the coordinate reference system.
#' @param treeBufferSize Size of the buffer added to each tree. Typically small but this seems to help when there
#'  is error in relative tree locations.
#' @param CHM Canopy height model in terra SpatRaster format.
#' @param searchRadius Maximum distance to offset plot location.
#' @param progress Boolean indicating whether or not progress dots should be printed.
#'
#' @return Data frame with search results.
#' @export
#'
#' @examples
testPlotLocations <- function(
  stemMap,
  coords = NULL,
  htField = NULL,
  plotX = 0,
  plotY = 0,
  plotRadius = 17.68,
  crs = NULL,
  treeBufferSize = 1.0,
  CHM = NULL,
  searchRadius = NULL,
  progress = TRUE
){
  # check parameters
  if (!isOpt(coords) || !isOpt(htField) || !isOpt(crs) || !isOpt(CHM) || !isOpt(searchRadius))
    stop("Invalid parameters")

  # build a grid of test locations
  offsetX <- seq(-searchRadius, searchRadius, by = terra::res(CHM)[1])
  offsetY <- seq(-searchRadius, searchRadius, by = terra::res(CHM)[1])

  # build data frame for output results
  res = data.frame(matrix(nrow = 0, ncol = 4))
  colnames(res) = c("offsetX", "offsetY", "correlation", "htError")

  for (j in 1:length(offsetX)) {
    for (k in 1:length(offsetY)) {
      # shift location of trees
      t <- stemMap[, c(coords, htField)]
      t[, coords[1]] <- t[, coords[1]] + offsetX[j]
      t[, coords[2]] <- t[, coords[2]] + offsetY[k]

      r <- computeCorrelationAndHeightError(
        t,
        coords = coords,
        htField = htField,
        plotX = plotX + offsetX[j],
        plotY = plotY + offsetY[k],
        plotRadius = plotRadius,
        crs = crs,
        treeBufferSize = treeBufferSize,
        CHM = CHM,
        searchRadius = searchRadius
      )

      # save offset and correlation between CHM and stems
      res <- rbind(res, data.frame("offsetX" = offsetX[j],
                                   "offsetY" = offsetY[k],
                                   "correlation" = r$correlation,
                                   "htError" = r$htError
                                    )
                    )
      if (progress) cat(".")
    }
    if (progress) cat("\n")
  }

  # compute the combined score
  res$combined <- (1.0 - (res$htError / max(res$htError))) * res$correlation

  return(res)
}

#' findBestPlotLocation
#'
#' @param searchResults Data frame returned by \code{testPlotLocations}.
#' @param rule Rule to be used to determine the best plot location. Possible
#' values are "maxCorrelation", "minHeightError", and "combined".
#'
#' @return Index into the data frame containing the offset from the original plot
#' location and the correlation and height error "scores".
#' @export
#'
#' @examples
findBestPlotLocation <- function(
    searchResults,
    rule = "maxcorrelation"
) {
  # get row for highest correlation
  if (tolower(rule) == "maxcorrelation")
    return(invisible(which.max(searchResults$corr)))
  else if (tolower(rule) == "minheighterror")
   return(invisible(which.min(searchResults$htError)))
  else if (tolower(rule) == "combined")
    return(invisible(which.max(searchResults$combined)))
  else stop(paste("Invalid rule:", rule))
}

#' rasterizeSearchResults
#'
#' Rasterize the results of a brute force search for plot locations.
#'
#' @param searchResults Data frame returned by \code{testPlotLocations}.
#' @param value Value used to populate the raster. Possible values are "correlation",
#' "heightError", and "combined".
#'
#' @return SpatRaster object
#' @export
#'
#' @examples
rasterizeSearchResults <- function(
    searchResults,
    value = "correlation"
) {
  r <- rast(xmin = min(searchResults$offsetX),
             ymin = min(searchResults$offsetY),
             xmax = max(searchResults$offsetX),
             ymax = max(searchResults$offsetY),
             resolution = searchResults$offsetY[2] - searchResults$offsetY[1])

  if (tolower(value) == "correlation")
    rr <- rasterize(as.matrix(searchResults[, 1:2]), r, values = searchResults[, 3])
  else if (tolower(value) == "heighterror")
    rr <- rasterize(as.matrix(searchResults[, 1:2]), r, values = searchResults[, 4])
  else if (tolower(value) == "combined")
    rr <- rasterize(as.matrix(searchResults[, 1:2]), r, values = searchResults[, 5])
  else
    stop(paste("Invalid value:", value))

  return(invisible(rr))
}