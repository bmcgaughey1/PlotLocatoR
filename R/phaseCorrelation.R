# useful links for the phase correlation approach
#
# Code to do image correlation to help adjust plot locations
# look for Fourier Mellin transform
# https://stackoverflow.com/questions/57801071/get-rotational-shift-using-phase-correlation-and-log-polar-transform
# https://sthoduka.github.io/imreg_fmt/docs/phase-correlation/
# https://sthoduka.github.io/imreg_fmt/docs/overall-pipeline/
# https://stackoverflow.com/questions/30630632/performing-a-phase-correlation-with-fft-in-r
# https://en.wikipedia.org/wiki/Phase_correlation
# in python https://stackoverflow.com/questions/2831527/phase-correlation



#'
#' xcoor3d: perform image registration using cross correlation. This function
#' was taken from the imagefx package and modified to better handle NA values
#' in the image matrices.
#'
#' @param img1 matrix: image representing reference set of objects
#' @param img2  matrix: image representing movable set of objects
#' @param normalize boolean: if TRUE, normalize the images by dividing
#'  by the mean value of each image.
#' @param addNoise boolean: if TRUE, add random noise to \code{img1}. This helps
#'  with the matching process when you have areas of high density trees with heights
#'  similar to those on the plot.
#' @param noiseMultiplier numeric value form 0.0 to 1.0 that is multiplied by the minimum
#'  and maximum values in \code{img1} to set limits on the values used to add noise to
#'  \code{img1}.
#'
#' @return \code{list} containing the following items:
#'    max.shifts shifts to align img2 with img1
#'    max.corr correlation values
#'    max.count
#'    corr.mat
#'
# ' @examples
xcorr3d <- function(
    img1,
    img2,
    normalize = TRUE,
    addNoise = TRUE,
    noiseMultiplier = 1.0
) {
  ## normalize by subtracting the mean
  if (normalize) {
    img1 <- img1-mean(img1, na.rm = TRUE)
    img2 <- img2-mean(img2, na.rm = TRUE)
  }

  if (addNoise) {
    img1 <- img1 + runif(length(img1), min(img1) * noiseMultiplier, max(img1) * noiseMultiplier)
  }

  ## go to the frequency domain and take conjugate of first image
  IMG1c <- Conj(stats::fft(img1))
  IMG2 <- stats::fft(img2)

  ## calculate the cross power spectrum
  R <- (IMG1c * IMG2)

  ## back to the spatial domain
  r.shift <- stats::fft(R,inverse=TRUE)/length(R)

  ## rearrange so zero frequency is in middle
  r <- imagefx::fftshift(Re(r.shift))

  ## normalize the r matrix
  r.norm <- r/(length(img1))/(stats::sd(as.vector(img1), na.rm = TRUE)*stats::sd(as.vector(img2), na.rm = TRUE))

  ## what is the correlation value associated with best shift
  max.cor <- max(r, na.rm = TRUE)

  max.count <- sum(r == max.cor)

  ## normalize the correlation value
  max.cor.norm <- max.cor/(length(img1))/(stats::sd(as.vector(img1), na.rm = TRUE)*stats::sd(as.vector(img2), na.rm = TRUE))

  ## find where the zero frequency is
  if(nrow(r)%%2 == 1) { zero.freq.x = -1.5 }
  if(nrow(r)%%2 == 0) { zero.freq.x = -1   }
  if(ncol(r)%%2 == 1) { zero.freq.y = -1.5 }
  if(ncol(r)%%2 == 0) { zero.freq.y = -1   }


  ## find the maximum value relative to the middle
  max.inds.abs <- which(r==max.cor,arr.ind=TRUE)-(dim(r)/2)

  ## adjust the maximum indices according to the true zero frequency
  max.inds <- max.inds.abs + c(zero.freq.x,zero.freq.y)

  ## if the input images are made of 0s, the r matrix could be filled
  ## with NaN and there is no max
  if(length(max.inds)==0){max.inds=matrix(c(0,0),nrow=1)}

  ## sometimes there is more than one max (normally in blank images)
  ## so take only the first max ind
  max.inds <- max.inds[1,]

  ## create a list to hold the max correlation value, its indices after
  ## shifting according to the zero frequency, and the original correlation matrix
  return.list = list()
  return.list$max.shifts <- max.inds
  return.list$max.corr <- max.cor.norm
  return.list$max.count <- max.count
  return.list$corr.mat <- r.norm

  ## return the above list
  return(return.list)
}


#' Compute the best location for the stem map relative to the CHM using cross or phase correlation.
#'
#' \code{findBestPlotLocationCorrelation} identifies tree approximate objects (TAOs) in the \code{CHM},
#' buffers the TAOs using \code{CHMbuffer}, and creates a raster image of the objects. \code{method}
#' controls the way the CHM is used for the image. When \code{method="buffer"}, the circular area
#' associated with each TAO is assigned the same height as the TAO. When \code{method="mask"}, the \code{CHM}
#' is masked using the TAO circles. The \code{stemMap} is used to create a second raster image representing
#' the trees on a plot. Each tree location is buffered using \code{stemMapBuffer} and all pixels in
#' the buffer are assigned the tree height.
#'
#' The shift of the stem map image that best aligns with the CHM image is determined using cross correlation
#' or phase correlation (depending on \code{corrMethod}) and converted into a new coordinate for the reference
#' point. This new coordinate can be used to translate the trees in the \code{stemMap}.
#'
#' @param CHM SpatRaster representing the canopy height surface
#' @param stemMap sf point object with tree positions based on \code{initialX,initialY} location
#' @param initialX initial reference point location
#' @param initialY initial reference point location
#' @param searchRadius radius defining area for possible locations. \code{searchRadius} will be used
#'  to crop the \code{CHM} to limit the area for the search. the \code{searchRadius} is only used when
#'  \code{cropCHM = TRUE}.
#' @param stemLocationFields character list with the field names in \code{stemMap} containing the X
#'  and Y values for the tree locations.
#' @param stemHeightField character name of the field in \code{stemMap} containing the tree height.
#' @param CHMbuffer radius for each TAO detected in \code{CHM}
#' @param stemMapBuffer radius for each tree in the \code{stemMap}
#' @param cropCHM boolean to indicate whether or not the CHM should be cropped to the search area
#'  defined by the \code{(initialX,initialY)} and \code{searchRadius}.
#' @param method character string defining the method used to handle the \code{CHM}. Possible
#'  values are "buffer" and "mask".
#' @param corrMethod character string specifying the correlation method used. Options are \code{"cross"}
#'  or \code{"phase"}.
#' @param normalize boolean: if TRUE, normalize the \code{CHM} and \code{stemMap} images by dividing
#'  by the mean value of each image.
#' @param addNoise boolean: if TRUE, add random noise to the \code{CHM} image when \code{corrMethod = "cross"}.
#'  This seems to help the matching process when you have areas of high density trees with heights
#'  similar to the tree heights on the plot. No noise is added to the \code{stemMap} image.
#' @param noiseMultiplier numeric value form 0.0 to 1.0 that is multiplied by the minimum
#'  and maximum values in the \code{CHM} image to set limits on the values used to add noise to
#'  the \code{CHM} image when using \code{addNoise = TRUE}. When using \code{addSelectiveNoise != "none"},
#'  the noise added to both  \code{CHM} and \code{stemMap} ranges from 1.0 to the maximum value in
#'  the \code{CHM} * \code{noiseMultiplier}.
#' @param addSelectiveNoise character string specifying which of the inputs should have random noise added to their
#'  background values. Choices are: "none", "both", "chm", and "stemmap". If not "none", add random noise to the
#'  \code{CHM} and/or \code{stemMap} image but only where values in the image(s) used for correlation have a value of 0.
#'  Like \code{addNoise}, this seems to help the matching process when you have areas
#'  of high density trees with heights similar to those on the plot. The best results seem to be obtained by setting
#'  \code{addSelectiveNoise = "chm"}.
#' @param includeRasters boolean: if TRUE, include the \code{CHM} and \code{stemMap} rasters in the
#'  returned list. This is primarily for debugging.
#'
#' @return invisible list containing the offsets from the \code{(initialX,initialY)} position
#'  and the (X,Y) for the best plot location.
#' @export
#'
# ' @examples
findBestPlotLocationCorrelation <- function(
    CHM,
    stemMap,
    initialX,
    initialY,
    searchRadius,
    stemLocationFields = c("X", "Y"),
    stemHeightField = "Ht",
    CHMbuffer = 1.0,
    stemMapBuffer = 1.0,
    cropCHM = TRUE,
    method = "buffer",
    corrMethod = "cross",
    normalize = TRUE,
    addNoise = FALSE,
    noiseMultiplier = 0.1,
    addSelectiveNoise = "chm",
    includeRasters = FALSE
) {
  NAReplaceValue <- 0

  # check noise options
  if (tolower(addSelectiveNoise) != "none" & addNoise)
    stop("You cannot use addNoise and addSelectiveNoise together")

  if (tolower(corrMethod) != "cross" & addNoise)
    stop("You cannot use addNoise when performing phase correlation")

#  if (tolower(corrMethod) == "cross" & tolower(addSelectiveNoise) != "none")
#    stop("Using AddSelectiveNoise when performing cross correlation has no effect")

  # crop CHM...if requested
  if (cropCHM) {
    CHM <- terra::crop(CHM
                , terra::ext(initialX - searchRadius, initialX + searchRadius, initialY - searchRadius, initialY + searchRadius)
                , snap = "near")
  }

  # replace NA values
  CHM[is.na(CHM)] <- NAReplaceValue

  # find TAOs
  # generate local maxima
  TAO <- lidR::locate_trees(CHM, lidR::lmf(ws = 5))
  TAO$Ht <- TAO$Z

  TAOxy <- sf::st_coordinates(TAO)

  TAO$X <- TAOxy[,1]
  TAO$Y <- TAOxy[,2]

  TAOtrees <- sf::st_as_sf(TAO,
                           coords = c("X", "Y"),
                           remove = FALSE,
                           crs = terra::crs(CHM))

  # buffer trees extracted from the CHM
  TAOtreesBuf <- sf::st_buffer(TAOtrees, CHMbuffer)

  # rasterize the buffered TAOs
  r <- terra::rast(terra::vect(TAOtreesBuf)
                   , extent = terra::ext(CHM)
                   , resolution = terra::res(CHM)[1]
                   , crs = terra::crs(CHM)
  )

  TAOPM <- terra::rasterize(terra::vect(TAOtreesBuf)
                            , r
                            , field = "Ht"
                            , extent = terra::ext(CHM)
                            , resolution = terra::res(CHM)[1]
                            , background = 0
  )

  if (tolower(method) == "mask") {
    # create CHM masked using TAOs extracted from CHM...this will keep the top of the crowns
    referenceCHM <- terra::mask(CHM, TAOPM, maskvalues = 0)
    referenceCHM[is.na(referenceCHM)] <- NAReplaceValue
  } else {
    # keep original CHM values within TAO buffers
    referenceCHM <- TAOPM
  }

  # create stem map object
  trees <- sf::st_as_sf(stemMap
                        , coords = stemLocationFields
                        , remove = FALSE
                        , crs = terra::crs(stemMap))

  # buffer trees...entire buffer around tree will have the same height
  treesBuf <- sf::st_buffer(trees, stemMapBuffer)

  r <- terra::rast(terra::vect(treesBuf)
                   , extent = terra::ext(CHM)
                   , resolution = terra::res(CHM)[1]
                   , crs = terra::crs(CHM)
  )

  PM <- terra::rasterize(terra::vect(treesBuf)
                         , r
                         , field = stemHeightField
                         , extent = terra::ext(CHM)
                         , resolution = terra::res(CHM)[1]
                         , background = 0
  )

  # add random noise for background cells
  #if (tolower(corrMethod) == "phase" & tolower(addSelectiveNoise) != "none") {
  if (tolower(addSelectiveNoise) != "none") {
      r <- terra::rast(extent = terra::ext(CHM)
                     , resolution = terra::res(CHM)[1]
                     , crs = terra::crs(CHM)
    )

    if (tolower(addSelectiveNoise) == "chm" | tolower(addSelectiveNoise) == "both") {
      # noise values
      #values(r) <- runif(nrow(r) * ncol(r), 1, minmax(referenceCHM)[[2]])
      values(r) <- runif(nrow(r) * ncol(r), 1, minmax(referenceCHM)[[2]] * noiseMultiplier)

      referenceCHM <- ifel(referenceCHM == 0, r, referenceCHM)
    }

    if (tolower(addSelectiveNoise) == "stemmap" | tolower(addSelectiveNoise) == "both") {
      # new noise values
      #values(r) <- runif(nrow(r) * ncol(r), 1, minmax(referenceCHM)[[2]])
      values(r) <- runif(nrow(r) * ncol(r), 1, minmax(referenceCHM)[[2]] * noiseMultiplier)

      PM <- ifel(PM == 0, r, PM)
    }
  }

  # create matrix objects
  i1 <- terra::as.matrix(referenceCHM, wide = TRUE)
  i2 <- terra::as.matrix(PM, wide = TRUE)

  # write matrix objects
  #saveRDS(i1, "G:/R_Stuff/ONRCDroneLidar/CHM_matrix.rds")
  #saveRDS(i2, "G:/R_Stuff/ONRCDroneLidar/PM_matrix.rds")

  # do correlation to get offset for stem map
  if (tolower(corrMethod) == "cross") {
    shifts <- xcorr3d(i1,i2, normalize = normalize, addNoise = addNoise)
  } else {
    shifts <- imagefx::pcorr3d(i1, i2)
  }

  # check for failure to find a better location
  if (shifts$max.shifts[1] == 0 & shifts$max.shifts[2] == 0) {
    offsetX <- 0
    offsetY <- 0
  } else {
    # shift vector is rotated 90 degree clockwise
    offsetX <- -shifts$max.shifts[2] * terra::res(CHM)[1]
    offsetY <- shifts$max.shifts[1] * terra::res(CHM)[1]
  }

  if (includeRasters) {
    invisible(return(list(offsetX = offsetX,
                          offsetY = offsetY,
                          newX = initialX + offsetX,
                          newY = initialY + offsetY,
                          CHM = referenceCHM,
                          stemMap = PM)
                      )
              )
  } else {
    invisible(return(list(offsetX = offsetX,
                          offsetY = offsetY,
                          newX = initialX + offsetX,
                          newY = initialY + offsetY)
                     )
              )
  }
}
