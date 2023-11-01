# function based on Gatziolis 2006

#' Compute the best location for the stem map relative to the CHM using the method outlined in Gatziolis (2006).
#'
#' \code{findBestPlotLocationGatziolis2006} identifies tree approximate objects (TAOs) in the \code{CHM}
#' and uses them to compute a pattern statistic that is compared to the pattern statistic computed using
#' the plot stem map. The statistic uses the azimuth from a center reference point (\code{CHM} grid point or
#' plot center), distance from the center reference point and tree height relative to the average height
#' for the search area or plot. The \code{CHM} grid point that produces the smallest difference in test
#' statistics is assumed to be the correct location for the plot.
#'
#' This methods works best in forests dominated by conifers and where there is some variability in tree
#' height. The methods performs poorly in areas where trees are uniformly spaced and have little variability
#' in height. Because the method relies on tree segmentation, it does not do well in forests where there is
#' a significant number of hardwood trees.
#'
#' Because this method does an exhaustive test of all cells in the \code{CHM}, it can be quite slow. Especially
#' when the \code{CHM} is large. You should almost always have \code{cropCHM = TRUE} to limit the area of the
#' search and set \code{searchRadius} to a value representing the largest expected error in \code{initialX,initialY}.
#'
#' @param CHM SpatRaster representing the canopy height surface
#' @param stemMap sf point object with tree positions based on \code{initialX,initialY} location. Attributes
#'  must include the location of each tree and the tree height. Column labels should be provided in
#'  \code{stemLocationFields} and \code{stemHeightField}.
#' @param initialX initial reference point location
#' @param initialY initial reference point location
#' @param plotRadius numeric value for the radius used to select trees for
#'  use when computing the test statistic. Trees with a distance less than or equal to \code{plotRadius}
#'  from \code{(initialX, initialY)} will be used to compute the test statistic.
#' @param searchRadius radius defining area for possible locations. \code{searchRadius} will be used
#'  to crop the \code{CHM} to limit the area for the search. The \code{searchRadius} is only used when
#'  \code{cropCHM = TRUE}.
#' @param stemLocationFields character list with the field names in \code{stemMap} containing the X
#'  and Y values for the tree locations.
#' @param stemHeightField character name of the field in \code{stemMap} containing the tree height.
#' @param TAOWindowSize numeric value specifying the diameter of the search window used to detect
#'  local maxima in the \code{CHM}. This is the \code{ws} parameter used with \code{lidR::locate_trees}.
#'  The value should be in the same horizontal units as the \code{CHM}.
#' @param CHMbuffer radius for each TAO detected in \code{CHM}
#' @param stemMapBuffer radius for each tree in the \code{stemMap}
#' @param cropCHM boolean to indicate whether or not the CHM should be cropped to the search area
#'  defined by the \code{(initialX,initialY)}, \code{searchRadius}, and \code{plotRadius}.
#' @param returnRaster boolean indicating that the raster containing the world view statistic should be returned.
#' @return invisible list containing the offsets from the \code{(initialX,initialY)} position
#'  and the (X,Y) for the best plot location. if \code{returnRaster=TRUE}, list also includes a SpatRaster
#'  object containing the world view statistic for each \code{CHM} cell in the search area. Cells outside
#'  the search area will have a value equal to the maximum value of the world view statistic.
#' @export
#'
# ' @examples
findBestPlotLocationGatziolis2006 <- function(
    CHM,
    stemMap,
    initialX,
    initialY,
    plotRadius,
    searchRadius,
    stemLocationFields = c("X", "Y"),
    stemHeightField = "Ht",
    TAOWindowSize = 5,
    cropCHM = TRUE,
    returnRaster = FALSE
) {
  NAReplaceValue <- 0

  # crop CHM...if requested
  if (cropCHM) {
    CHM <- terra::crop(CHM
                       , terra::ext(initialX - searchRadius - plotRadius,
                                    initialX + searchRadius + plotRadius,
                                    initialY - searchRadius - plotRadius,
                                    initialY + searchRadius + plotRadius)
                       , snap = "near")
  }

  # replace NA values
  CHM[is.na(CHM)] <- NAReplaceValue

  # find TAOs
  # generate local maxima
  TAO <- lidR::locate_trees(CHM, lidR::lmf(ws = TAOWindowSize))
  TAO$Ht <- TAO$Z

  TAOxy <- sf::st_coordinates(TAO)

  TAO$X <- TAOxy[,1]
  TAO$Y <- TAOxy[,2]

  TAO <- sf::st_drop_geometry(TAO)

  # compute test statistic for plot
  ps <- computeWorldViewStatistic(stemMap, initialX, initialY, plotRadius, stemLocationFields, stemHeightField)

  # create raster to store results
  r <- terra::rast(extent = terra::ext(CHM)
                   , resolution = terra::res(CHM)[1]
                   , crs = terra::crs(CHM)
                   , vals = 0
  )

  # walk through cells in CHM and compute test statistic for each grid point and compare to plot
  # ***** might be able to do this with app() but not sure how to get the xy
  for (row in 1:nrow(CHM)) {
    if ((row %% 5) == 0) message("row:", row, " of ", nrow(CHM))
    for (col in 1:ncol(CHM)) {

      X <- terra::xyFromCell(CHM, terra::cellFromRowCol(CHM, row, col))[[1]]
      Y <- terra::xyFromCell(CHM, terra::cellFromRowCol(CHM, row, col))[[2]]

      if (sqrt((initialX - X) * (initialX - X) + (initialY - Y) * (initialY - Y)) <= searchRadius) {
        s <- computeWorldViewStatistic(TAO, X, Y, plotRadius, c("X", "Y"), "Ht")

        diff <- ps - s

        r[row, col] <- r[row, col] + sum(abs(diff))
      }
    }
  }

  # set cells with -1 to large value
  r <- terra::classify(r, cbind(0, terra::minmax(r)[[2]]))

  # find CHM cell with minimum difference...could be multiple cells
  c <- terra::where.min(r, list = TRUE)

  # get the first cell number with the minimum value
  c <- c[[1]][[1]]

  newX <- terra::xyFromCell(r, c)[1]
  newY <- terra::xyFromCell(r, c)[2]
  offsetX <- newX - initialX
  offsetY <- newY - initialY

  invisible(
    if (returnRaster) {
      return(list(offsetX = offsetX,
                  offsetY = offsetY,
                  newX = initialX + offsetX,
                  newY = initialY + offsetY,
                  wv = r)
      )
    }
    else {
      return(list(offsetX = offsetX,
                  offsetY = offsetY,
                  newX = initialX + offsetX,
                  newY = initialY + offsetY)
      )
    }
  )
}

#' Compute the test statistic for a stem map using a center point and radius
#'
#' \code{computeWorldViewStatistic} can be used for plot stem maps or a larger stem map.
#'
#' @param stemMap sf point object with tree positions based on \code{initialX,initialY} location. Attributes
#'  must include the location of each tree and the tree height. Column labels should be provided in
#'  \code{stemLocationFields} and \code{stemHeightField}.
#' @param centerX numeric value for the reference point X
#' @param centerY numeric value for the reference point Y
#' @param plotRadius numeric value for the radius used to select trees from the \code{stemMap} for
#'  use when computing the test statistic. Trees with a distance less than or equal to \code{plotRadius}
#'  from \code{(centerX, centerY)} will be used to compute the test statistic.
#' @param stemLocationFields character list with the field names in \code{stemMap} containing the X
#'  and Y values for the tree locations.
#' @param stemHeightField character name of the field in \code{stemMap} containing the tree height.
#'
#' @return invisible numeric value for the world view test statistic.
#' @export
#'
#' @examples
computeWorldViewStatistic <- function(
  stemMap,
  centerX,
  centerY,
  plotRadius,
  stemLocationFields = c("X", "Y"),
  stemHeightField = "Ht"
) {
  # compute distance from (centerX,centerY) to each tree
  stemMap$ttt_distance <- sqrt((stemMap[, stemLocationFields[1]] - centerX) * (stemMap[, stemLocationFields[1]] - centerX) + (stemMap[, stemLocationFields[2]] - centerY) * (stemMap[, stemLocationFields[2]] - centerY))

  # filter to get trees within plotRadius
  stemMap <- stemMap[stemMap$ttt_distance <= plotRadius, ]

  # compute azimuth for test trees
  stemMap$ttt_angle <- ((atan2(stemMap[, stemLocationFields[2]] - centerY, stemMap[, stemLocationFields[1]] - centerX) * (180.0 / pi)) + 360) %% 360

#  print(stemMap)

  # compute mean height and variance of heights
  mean <- mean(stemMap[, stemHeightField])
  var <- var(stemMap[, stemHeightField])

  # initialize vector of test values
  ts <- vector("numeric", 360)

  # compute vector of test values
  for (i in 1:nrow(stemMap)) {
    ts[floor(stemMap$ttt_angle[i])] <- ts[floor(stemMap$ttt_angle[i])] + (plotRadius - stemMap$ttt_distance[i] + abs((stemMap[i, stemHeightField] - mean) / mean) * var)
  }

  # return test values
  invisible(return(ts))
}
